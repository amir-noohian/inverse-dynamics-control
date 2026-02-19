/*
 * 4dof_gravcomp_zeus.cpp 
 *
 *  Created on: July, 2024
 *      Author: Amir Noohian
 */


#include <dynamics_4dof_slax.hpp>
#include <constant_vel_refference_traj.hpp>
#include <sin_jp_trajectory.hpp>
#include <unistd.h>
#include <iostream>
#include <string>
#include <barrett/units.h>
#include <barrett/systems.h>
#include <barrett/products/product_manager.h>
#include <barrett/detail/stl_utils.h>
#include <barrett/log.h>
#include <barrett/standard_main_function.h>

using namespace barrett;

//The WAM arm has a torque saturation limit
template <size_t DOF>
typename units::JointTorques<DOF>::type saturateJt(
    const typename units::JointTorques<DOF>::type& x,
    const typename units::JointTorques<DOF>::type& limit) {
    int index;
    double minRatio;

	minRatio = limit.cwiseQuotient(x.cwiseAbs()).minCoeff(&index);
    if (minRatio < 1.0) {
        return minRatio * x;
    } else {
        return x;
    }
}

//Our joint space inverse dynamics control (idc) : id_feedFWD + PD controller
template<size_t DOF>
class JsIDController :  public systems::System{ //
	BARRETT_UNITS_TEMPLATE_TYPEDEFS(DOF);

public:
	Input<jt_type> gravityInput; //Gravity
	Input<jt_type> feedFWDInput; //Inverse Dynamic Feed FWD

// IO  (outputs)
public:
	Output<jt_type> controlJtOutput;    // Output Control Joint Torque

protected:
	typename Output<jt_type>::Value* controlJtOutputValue;

public:
jt_type computedT;

public:
	explicit JsIDController(const std::string& sysName = "JsIDController"):
		System(sysName), controlJtOutput(this, &controlJtOutputValue), feedFWDInput(this), gravityInput(this)/*, kpCf(kpCoeff), kdCf(KdCoeff)*/{}

	virtual ~JsIDController() { this->mandatoryCleanUp(); }

protected:
	jt_type jt_out, feedFWD, gravity, gravity_mod;

	virtual void operate() {

		/*Taking feed forward term from the input terminal of this system*/
		gravity = this->gravityInput.getValue();
		feedFWD = this->feedFWDInput.getValue();

		gravity_mod = gravity;
        gravity_mod[4] = 0.0;
        gravity_mod[5] = 0.0;
        gravity_mod[6] = 0.0;

		jt_out = feedFWD - gravity_mod;

		computedT = jt_out;
		
		this->controlJtOutputValue->setData(&jt_out);
	}

private:
	DISALLOW_COPY_AND_ASSIGN(JsIDController);
};

//Creating a templated multiplier for our real-time computation
template<typename T1, typename T2, typename OutputType>
  class Multiplier : public systems::System, public systems::SingleOutput<OutputType>
  {
  public:
    Input<T1> input1;
  public:
    Input<T2> input2;

  OutputType Data;
  public:
    Multiplier(std::string sysName = "Multiplier") :
        systems::System(sysName), systems::SingleOutput<OutputType>(this), input1(this), input2(this)
    {
    }
    virtual ~Multiplier()
    {
      mandatoryCleanUp();
    }

  protected:
    OutputType data;
    virtual void operate()
    {
      data = input1.getValue() * input2.getValue();
      Data = data;
      this->outputValue->setData(&data);
    }

  private:
    DISALLOW_COPY_AND_ASSIGN(Multiplier);

  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  };

template<size_t DOF>
int wam_main(int argc, char** argv, ProductManager& pm,	systems::Wam<DOF>& wam) {
    BARRETT_UNITS_TEMPLATE_TYPEDEFS(DOF);

	char tmpFile[] = "/tmp/btXXXXXX";
	if (mkstemp(tmpFile) == -1) {
		printf("ERROR: Couldn't create temporary file!\n");
		return 1;
	}

	pm.getSafetyModule()->setVelocityLimit(1.5);


	wam.gravityCompensate();

	JsIDController<DOF> IDController;
	Dynamics<DOF> wam4dofDynamics;
	
	jt_type jtLimits;
	jtLimits << 50, 40, 30, 30, 10, 10, 10;
	systems::Callback<jt_type> jtSat(boost::bind(saturateJt<DOF>,_1, 0.5*jtLimits));

	systems::Ramp time(pm.getExecutionManager(), 1.0);


	jv_type jv;
    jv.setConstant(0.0);
    systems::Constant<jv_type> zeroVelocity(jv);
    pm.getExecutionManager()->startManaging(zeroVelocity);

	ja_type ja;
    ja.setConstant(0.0);
    systems::Constant<ja_type> zeroAcceleration(ja);
    pm.getExecutionManager()->startManaging(zeroAcceleration);

	systems::connect(wam.jpOutput, wam4dofDynamics.jpInputDynamics);
	systems::connect(zeroVelocity.output, wam4dofDynamics.jvInputDynamics);
	systems::connect(zeroAcceleration.output, wam4dofDynamics.jaInputDynamics);
	systems::connect(wam4dofDynamics.dynamicsFeedFWD, IDController.feedFWDInput);
	systems::connect(wam.gravity.output, IDController.gravityInput);
	systems::connect(IDController.controlJtOutput, jtSat.input);


//	RT Logging stuff
	systems::Ramp timelog(pm.getExecutionManager(), 1.0);
	systems::TupleGrouper<double, jp_type, jt_type, jt_type, jt_type> tg;
	systems::connect(timelog.output, tg.template getInput<0>());
	systems::connect(wam.jpOutput, tg.template getInput<1>());
	systems::connect(wam.jtSum.output, tg.template getInput<2>());
	systems::connect(wam.gravity.output, tg.template getInput<3>());
	systems::connect(wam4dofDynamics.dynamicsFeedFWD, tg.template getInput<4>());


	typedef boost::tuple<double, jp_type, jt_type, jt_type, jt_type> tuple_type;
	const size_t PERIOD_MULTIPLIER = 1;
	systems::PeriodicDataLogger<tuple_type> logger(
			pm.getExecutionManager(),
			new log::RealTimeWriter<tuple_type>(tmpFile, PERIOD_MULTIPLIER * pm.getExecutionManager()->getPeriod()),
			PERIOD_MULTIPLIER);
	
	time.stop();
	time.reset();

	timelog.stop();
	timelog.reset();
	timelog.start();
	systems::connect(tg.output, logger.input);
	printf("Logging started.\n");
	
	std::cout<<"Press [Enter] to start."<<std::endl;
	detail::waitForEnter();
	time.start();
	usleep(2500);
	wam.trackReferenceSignal(jtSat.output); // Put it after time start so that feefwd term is active when we start applying torque.


	std::cout<<"Press [Enter] to stop."<<std::endl;
	detail::waitForEnter();
	time.stop();

	timelog.stop();
	logger.closeLog();
	printf("Logging stopped.\n");

	log::Reader<tuple_type> lr(tmpFile);
	lr.exportCSV(argv[1]);
	printf("Output written to %s.\n", argv[1]);
	std::remove(tmpFile);

	wam.moveHome();

	// Wait for the user to press Shift-idle
	pm.getSafetyModule()->waitForMode(SafetyModule::IDLE);
	return 0;
}