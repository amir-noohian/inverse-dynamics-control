/*
 * 2DOF_wam_joint_space_ctc.cpp 
 *
 *  Created on: Aug, 2024
 *      Author: Faezeh
 */


#include <Gravity.hpp>
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

    minRatio = (limit.cwise() / (x.cwise().abs())).minCoeff(&index);
    if (minRatio < 1.0) {
        return minRatio * x;
    } else {
        return x;
    }
}

//Our joint space Computed Torque Controller(CTC) : id_feedFWD + gravity + PD controller
template<size_t DOF>
class JsCTCController :  public systems::System{ //
	BARRETT_UNITS_TEMPLATE_TYPEDEFS(DOF);

// IO  (inputs)
public:
	Input<jt_type> feedFWDInput; //Inverse Dynamic Feed FWD
	Input<jt_type> gravityInput; //Gravity Term Input from WAM

// IO  (outputs)
public:
	Output<jt_type> controlJtOutput;    // Output Control Joint Torque

protected:
	typename Output<jt_type>::Value* controlJtOutputValue;

public:
jt_type computedT;

public:
	explicit JsCTCController(const std::string& sysName = "JsCTCController"):
		System(sysName), controlJtOutput(this, &controlJtOutputValue), feedFWDInput(this), gravityInput(this) /*, kpCf(kpCoeff), kdCf(KdCoeff)*/{}

	virtual ~JsCTCController() { this->mandatoryCleanUp(); }

protected:
	jt_type jt_out, feedFWD, gravity;

	virtual void operate() {
		/*Taking feed forward term from the input terminal of this system*/
		feedFWD = this->feedFWDInput.getValue();
		gravity = this->gravityInput.getValue();

		gravity[1] = 0.0;
		gravity[3] = 0.0;

		jt_out = feedFWD + gravity;

		computedT = jt_out;
		
		this->controlJtOutputValue->setData(&jt_out);
	}

private:
	DISALLOW_COPY_AND_ASSIGN(JsCTCController);
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

	JsCTCController<DOF> compTorqueController;//, inputs1[0], inputs1[1]);
	Dynamics<DOF> wam2dofDynamics;
	jt_type jtLimits;
	jtLimits << 25, 20, 15, 15, 5, 5, 5;
	systems::Callback<jt_type> jtSat(boost::bind(saturateJt<DOF>,_1, 2*jtLimits));
	systems::Ramp time(pm.getExecutionManager(), 1.0);
    
    systems::connect(wam.jpController, wam2dofDynamics.jpInputDynamics);
   	systems::connect(wam.gravity.output, compTorqueController.gravityInput);
    systems::connect(wam2dofDynamics.dynamicsFeedFWD, compTorqueController.feedFWDInput);
    systems::connect(compTorqueController.controlJtOutput, jtSat.input);

//	RT Logging stuff
	systems::Ramp timelog(pm.getExecutionManager(), 1.0);
	systems::TupleGrouper<double, jp_type, jt_type, jt_type, jt_type, jt_type> tg;
	systems::connect(timelog.output, tg.template getInput<0>());
	systems::connect(wam.jpOutput, tg.template getInput<1>());
	systems::connect(wam.jtSum.output, tg.template getInput<2>());
	systems::connect(wam.gravity.output, tg.template getInput<3>());
	systems::connect(compTorqueController.controlJtOutput, tg.template getInput<4>());
	systems::connect(wam2dofDynamics.dynamicsFeedFWD, tg.template getInput<5>());
	// systems::connect(refTraj.jaDes, tg.template getInput<7>());

	typedef boost::tuple<double, jp_type, jt_type, jt_type, jt_type, jt_type> tuple_type;
	const size_t PERIOD_MULTIPLIER = 1;
	systems::PeriodicDataLogger<tuple_type> logger(
			pm.getExecutionManager(),
			new log::RealTimeWriter<tuple_type>(tmpFile, PERIOD_MULTIPLIER * pm.getExecutionManager()->getPeriod()),
			PERIOD_MULTIPLIER);
	
	timelog.stop();
	timelog.reset();
	timelog.start();
	systems::connect(tg.output, logger.input);
	printf("Logging started.\n");
	
	wam.trackReferenceSignal(jtSat.output); // Put it after time start so that feefwd term is active when we start applying torque.

	// Wait for the user to press Shift-idle
	pm.getSafetyModule()->waitForMode(SafetyModule::IDLE);

	timelog.stop();
	logger.closeLog();
	printf("Logging stopped.\n");

	log::Reader<tuple_type> lr(tmpFile);
	lr.exportCSV(argv[1]);
	printf("Output written to %s.\n", argv[1]);
	std::remove(tmpFile);
	return 0;
}