/*
 * 4dof_id_control.cpp 
 *
 *  Created on: Juen, 2024
 *      Author: Amir Noohian
 */


#include <dynamics_4dof_old.hpp>
#include <dynamics_4dof.hpp>
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

// IO  (inputs)
public:
	Input<jp_type> refJPInput;	// reference joint pos. input   
	Input<jv_type> refJVInput;	// reference joint vel. input

public:
	Input<jp_type> feedbackjpInput;     // joint pos. input       
	Input<jv_type> feedbackjvInput;  	// joint vel. input

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
	explicit JsIDController(Eigen::MatrixXd proportionalGains, Eigen::MatrixXd dampingGains/*, float kpCoeff, float KdCoeff*/, const std::string& sysName = "JsIDController"):
		System(sysName), refJPInput(this), refJVInput(this), feedbackjpInput(this), feedbackjvInput(this), controlJtOutput(this, &controlJtOutputValue), feedFWDInput(this), gravityInput(this), kp(proportionalGains), kd(dampingGains)/*, kpCf(kpCoeff), kdCf(KdCoeff)*/{}

	virtual ~JsIDController() { this->mandatoryCleanUp(); }

protected:
	float kpCf, kdCf; 
	sqm_type kp, kd;
	jt_type jt_out, feedFWD, gravity;
	jp_type jp_sys, jp_ref;
	jv_type jv_sys, jv_ref;
	ja_type ja_ref;

	virtual void operate() {
		/*Taking reference values from the input terminal of this system*/
		jp_ref = this->refJPInput.getValue();
		jv_ref = this->refJVInput.getValue();

		/*Taking feedback values from the input terminal of this system*/
		jp_sys = this->feedbackjpInput.getValue();
		jv_sys = this->feedbackjvInput.getValue();

		/*Taking feed forward term from the input terminal of this system*/
		feedFWD = this->feedFWDInput.getValue();
		gravity = this->gravityInput.getValue();

		jt_out = 1 * feedFWD - 1 * gravity +  kp * (jp_ref - jp_sys) + kd * (jv_ref - jv_sys);

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

	//Moving to start pose
	jp_type start_pos;
	start_pos<< 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;

	wam.gravityCompensate();
	usleep(250000);
	wam.moveTo(start_pos);

	//Trapezoidal vel traj
	// jp_type end_pos;
	// end_pos<< 0.0, 0.0, 0.0, -0.5, 0.0, 0.0, 0.0;
	// std::vector<jp_type, Eigen::aligned_allocator<jp_type> > vec;
	// vec.push_back(start_pos);
	// vec.push_back(end_pos);
	// math::Spline<jp_type> spline(vec);
	// math::TrapezoidalVelocityProfile profile(0.0, 0.0, 0.0, spline.changeInS());
	// systems::Callback<double, jp_type> refTraj(boost::bind(boost::ref(spline), boost::bind(boost::ref(profile), _1)));

	// Sinosoidal Traj
	// std::vector<float> inputs;
	// float temp;
	// std::cout<<"Enter A2, A4, and f"<<std::endl;
	
	// for (int i = 0; i < 3; ++i) {
    //     std::cin >> temp;
    //     inputs.push_back(temp);  // Add the number to the vector
    // }

	v_type A;
	A << 0.0, 0.4 , 0.0, -0.4, 0.0, 0.0, 0.0;
	double f = 0.1;
	sinJpRefTrajectory<DOF> refTraj(start_pos, A, f);

	// //Const Vel Trajectory
	// jv_type thetad_des;
    // thetad_des << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
	// std::vector<float> inputs;
	// float temp;
	// std::cout<<"Enter setpoit 2 and 4"<<std::endl;
	// for (int i = 0; i < 2; ++i) {
    //     std::cin >> temp;
    //     inputs.push_back(temp);  // Add the number to the vector
    // }
	// start_pos << 0.0, inputs[0], 0.0, inputs[1], 0.0, 0.0, 0.0;
    // constVelRefTrajectory<DOF> refTraj(start_pos, thetad_des);

	// float kpCf, kdCf;
	v_type kp,kd; 
	// std::vector<float> inputs1;
	// std::cout<<"Enter kp 2 and 4, kd 2 and 4"<<std::endl;
	// for (int i = 0; i < 4; ++i) {
    //     std::cin >> temp;
    //     inputs1.push_back(temp);  // Add the number to the vector
    // }
	
	kp << 900, 2500,  600,  500,   50,   50,    8;
	kd << 10,   20,    5,    2,  0.5,  0.5, 0.05;
	// 	kp << 900, 2500,  600,  500,   50,   50,    8; //250 50
	// kd << 10,   20,    5,    2,  0.5,  0.5, 0.05; //10 1
	JsIDController<DOF> IDController(kp.asDiagonal(),kd.asDiagonal());//, inputs1[0], inputs1[1]);
	Dynamics<DOF> wam4dofDynamics;
	// DynamicsOld<DOF> wam4dofDynamicsOld;
	jt_type jtLimits;
	jtLimits << 25, 20, 15, 15, 5, 5, 5;
	systems::Callback<jt_type> jtSat(boost::bind(saturateJt<DOF>,_1, 1*jtLimits));
	systems::Ramp time(pm.getExecutionManager(), 1.0);

	systems::PrintToStream<jt_type> printID(pm.getExecutionManager(), "ID: ");

    
	// double omega_p = 1000.0;
	// systems::FirstOrderFilter<jp_type> hp1;
	// hp1.setHighPass(jp_type(omega_p), jp_type(omega_p));
	// systems::FirstOrderFilter<jp_type> hp2;
	// hp2.setHighPass(jp_type(omega_p), jp_type(omega_p));
	// systems::Gain<jp_type, double, jv_type> jvDes(1.0);
	// systems::Gain<jp_type, double, ja_type> jaDes(1.0);

	systems::connect(time.output, refTraj.timef);
	//systems::connect(time.output, refTraj.input);
	// connect(trefTraj.output, hp1.input);
	// connect(hp1.output, hp2.input);
	// connect(hp2.output, jaDes.input);
	// connect(hp1.output, jvDes.input);
	// pm.getExecutionManager()->startManaging(hp2);
	// sleep(1);

	// systems::connect(refTraj.output, compTorqueController.refJPInput);
    // systems::connect(jvDes.output, compTorqueController.refJVInput);	

    // systems::connect(refTraj.output, wam2dofDynamics.jpInputDynamics);
	// systems::connect(jvDes.output, wam2dofDynamics.jvInputDynamics);
    // systems::connect(jaDes.output, wam2dofDynamics.jaInputDynamics);

    systems::connect(refTraj.referencePTrack, IDController.refJPInput);
    systems::connect(refTraj.referenceVTrack, IDController.refJVInput);	

    systems::connect(refTraj.referencePTrack, wam4dofDynamics.jpInputDynamics);
	systems::connect(refTraj.referenceVTrack, wam4dofDynamics.jvInputDynamics);
    systems::connect(refTraj.referenceATrack, wam4dofDynamics.jaInputDynamics);

    systems::connect(wam4dofDynamics.dynamicsFeedFWD, IDController.feedFWDInput);
	systems::connect(wam.gravity.output, IDController.gravityInput);
	systems::connect(wam.jpOutput, IDController.feedbackjpInput);
	systems::connect(wam.jvOutput, IDController.feedbackjvInput);
    systems::connect(IDController.controlJtOutput, jtSat.input);

	// systems::connect(wam4dofDynamics.dynamicsFeedFWD, printID.input);

	// systems::connect(refTraj.referencePTrack, wam4dofDynamics.jpInputDynamics);
	// systems::connect(refTraj.referenceVTrack, wam4dofDynamics.jvInputDynamics);
    // systems::connect(refTraj.referenceATrack, wam4dofDynamics.jaInputDynamics);

//	RT Logging stuff
	systems::Ramp timelog(pm.getExecutionManager(), 1.0);
	systems::TupleGrouper<double, jp_type, jp_type, jt_type, jt_type, jt_type, jt_type, jv_type, ja_type> tg;
	systems::connect(timelog.output, tg.template getInput<0>());
	systems::connect(refTraj.referencePTrack, tg.template getInput<1>());
	systems::connect(wam.jpOutput, tg.template getInput<2>());
	systems::connect(wam.jtSum.output, tg.template getInput<3>());
	systems::connect(wam.gravity.output, tg.template getInput<4>());
	systems::connect(IDController.controlJtOutput, tg.template getInput<5>());
	systems::connect(wam4dofDynamics.dynamicsFeedFWD, tg.template getInput<6>());
	systems::connect(refTraj.referenceVTrack, tg.template getInput<7>());
	systems::connect(refTraj.referenceATrack, tg.template getInput<8>());
	// systems::connect(refTraj.jaDes, tg.template getInput<7>());

	typedef boost::tuple<double, jp_type, jp_type, jt_type, jt_type, jt_type, jt_type, jv_type, ja_type> tuple_type;
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
	// wam.gravityCompensate(false);
	time.start();
	usleep(2500);
	wam.trackReferenceSignal(jtSat.output); // Put it after time start so that feefwd term is active when we start applying torque.

	// while (pm.getSafetyModule()->getMode() == SafetyModule::ACTIVE){
	// 	// std::cout <<"Control output:"<< compTorqueController.computedT << std::endl;
	// 	// if(compTorqueController.feedFWDInput.valueDefined()) std::cout<<"Feedforward input:"<<compTorqueController.feedFWDInput.getValue()<<std::endl;
	// 	if(compTorqueController.gravityInput.valueDefined()) std::cout<<"Gravity:"<<compTorqueController.gravityInput.getValue()<<std::endl;
	// 	if(jtSat.input.valueDefined()) std::cout<<"Control output:"<<jtSat.input.getValue()<<std::endl;
	// 	// if(trajectory.input.valueDefined()) std::cout<<"Time:"<<trajectory.input.getValue()<<std::endl;
	// 	if(compTorqueController.refJPInput.valueDefined()) std::cout<<"Reference position:"<<compTorqueController.refJPInput.getValue()<<std::endl;
	// 	if(compTorqueController.feedFWDInput.valueDefined()) std::cout<<"Feed FWD:"<<compTorqueController.feedFWDInput.getValue()<<std::endl;

	// 	std::cout<<"Pos:"<<wam.getJointPositions()<<std::endl;
	// 	std::cout<<"Vel:"<<wam.getJointVelocities()<<std::endl;

	// 	//btsleep(0.2);
	// }

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