/*
 * inverse_dynamics_test.cpp
 *
 *  Created on: Nov 8, 2011
 *      Author: dc
 */
#include <Dynamics.hpp>
#include <cstdlib>  // For mkstmp()
#include <cstdio>  // For remove()

#include <boost/tuple/tuple.hpp>

#define EIGEN_USE_NEW_STDVECTOR
#include <Eigen/StdVector>

#include <libconfig.h++>

#include <barrett/log.h>
#include <barrett/units.h>
#include <barrett/systems.h>
#include <barrett/products/product_manager.h>

#define BARRETT_SMF_VALIDATE_ARGS
#include <barrett/standard_main_function.h>


using namespace barrett;
using systems::connect;


bool validate_args(int argc, char** argv) {
	if (argc != 2) {
		printf("Usage: %s <fileName>\n", argv[0]);
		return false;
	}
	return true;
}

template<size_t DOF>
int wam_main(int argc, char** argv, ProductManager& pm, systems::Wam<DOF>& wam) {
	BARRETT_UNITS_TEMPLATE_TYPEDEFS(DOF);

	char tmpFile[] = "/tmp/btXXXXXX";
	if (mkstemp(tmpFile) == -1) {
		printf("ERROR: Couldn't create temporary file!\n");
		return 1;
	}

	//Moving to start pose
	jp_type start_pos;
	start_pos<< 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;

	wam.gravityCompensate();
	usleep(250000);

	wam.moveTo(start_pos);

	//Trapezoidal vel traj
	jp_type end_pos;
	end_pos<< 0.0, 0.0, 0.0, -0.5, 0.0, 0.0, 0.0;
	std::vector<jp_type, Eigen::aligned_allocator<jp_type> > vec;
	vec.push_back(start_pos);
	vec.push_back(end_pos);
	math::Spline<jp_type> spline(vec);
	math::TrapezoidalVelocityProfile profile(0.5, 0.5, 0.0, spline.changeInS());
	systems::Callback<double, jp_type> refTraj(boost::bind(boost::ref(spline), boost::bind(boost::ref(profile), _1)));
	systems::Ramp moveTime(pm.getExecutionManager(), 1.0);

	connect(moveTime.output, refTraj.input);
	wam.trackReferenceSignal(refTraj.output);

	double omega_p = 1000.0;
	systems::FirstOrderFilter<jp_type> hp1;
	hp1.setHighPass(jp_type(omega_p), jp_type(omega_p));
	systems::FirstOrderFilter<jp_type> hp2;
	hp2.setHighPass(jp_type(omega_p), jp_type(omega_p));
	systems::Gain<jp_type, double, jv_type> jvDes(1.0);
	systems::Gain<jp_type, double, ja_type> jaDes(1.0);
    Dynamics<DOF> feedFWD;

	connect(refTraj.output, hp1.input);
	connect(hp1.output, hp2.input);
	connect(hp2.output, jaDes.input);
	connect(hp1.output, jvDes.input);
	pm.getExecutionManager()->startManaging(hp2);
	sleep(1);

	connect(refTraj.output, feedFWD.jpInputDynamics);
	connect(jvDes.output, feedFWD.jvInputDynamics);
    connect(jaDes.output, feedFWD.jaInputDynamics);
	connect(feedFWD.dynamicsFeedFWD, wam.input);

	wam.jpController.setKp(wam.jpController.getKp() * 0.1);
	// wam.jpController.setKp(wam.jpController.getKp() * 0.1);
	wam.jpController.setKd(wam.jpController.getKd() * 0.5);


// RT Logging stuff
	systems::Ramp timelog(pm.getExecutionManager(), 1.0);
	systems::TupleGrouper<double, jp_type, jp_type, jt_type, jt_type, jt_type, jt_type, jv_type, ja_type> tg;
	systems::connect(timelog.output, tg.template getInput<0>());
	systems::connect(refTraj.output, tg.template getInput<1>());
	systems::connect(wam.jpOutput, tg.template getInput<2>());
	systems::connect(wam.jtSum.output, tg.template getInput<3>());
	systems::connect(wam.gravity.output, tg.template getInput<4>());
	systems::connect(wam.gravity.output, tg.template getInput<5>());
	systems::connect(feedFWD.dynamicsFeedFWD, tg.template getInput<6>());
	systems::connect(jvDes.output, tg.template getInput<7>());
	systems::connect(jaDes.output, tg.template getInput<8>());

	typedef boost::tuple<double, jp_type, jp_type, jt_type, jt_type, jt_type, jt_type, jv_type, ja_type> tuple_type;
	const size_t PERIOD_MULTIPLIER = 1;
	systems::PeriodicDataLogger<tuple_type> logger(
			pm.getExecutionManager(),
			new log::RealTimeWriter<tuple_type>(tmpFile, PERIOD_MULTIPLIER * pm.getExecutionManager()->getPeriod()),
			PERIOD_MULTIPLIER);

	timelog.start();
	connect(tg.output, logger.input);
	printf("Logging started.\n");

	usleep(100000);
	moveTime.start();
	sleep(4);

	logger.closeLog();
	printf("Logging stopped.\n");

	log::Reader<tuple_type> lr(tmpFile);
	lr.exportCSV(argv[1]);
	printf("Output written to %s.\n", argv[1]);
	std::remove(tmpFile);

	wam.moveHome();
	pm.getSafetyModule()->waitForMode(SafetyModule::IDLE);

	return 0;
}
