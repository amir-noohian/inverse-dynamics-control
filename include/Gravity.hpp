#ifndef DYNAMICS_H_
#define DYNAMICS_H_

#include <barrett/math/traits.h>
#include <list>
#include <barrett/units.h>
#include <barrett/detail/ca_macro.h>
#include <barrett/systems/abstract/system.h>

#include <eigen3/Eigen/Core>
#include <libconfig.h++>

#include <barrett/detail/ca_macro.h>
#include <barrett/math/traits.h>
#include <barrett/systems/abstract/execution_manager.h>
#include <barrett/systems/abstract/controller.h>

#include <PI_2D_onlygravity.hpp>
#include <Y_2D_onlygravity.hpp>

using namespace barrett;

template<size_t DOF>
class Dynamics: public systems::System {

	BARRETT_UNITS_TEMPLATE_TYPEDEFS(DOF);

	/* Torque*/
public:
	Input<jp_type> jpInputDynamics;

public:
	Output<jt_type> dynamicsFeedFWD;
protected:
	typename Output<jt_type>::Value* dynamicsFeedFWDValue;

public:
	explicit Dynamics(/*systems::ExecutionManager* em*/) :
			jpInputDynamics(this), dynamicsFeedFWD(this, &dynamicsFeedFWDValue) {
	//		      em->startManaging(*this);
//		    }
	}

	virtual ~Dynamics() {
		this->mandatoryCleanUp();
	}

protected:
	jp_type tmp_theta_pos;
	jv_type tmp_theta_vel;
	ja_type tmp_theta_acc;
	Eigen::Matrix<double,2 ,4> Y;
	Eigen::Matrix<double, 4, 1> P;
	Eigen::Vector4d ThetaInput;
	Eigen::Vector4d ThetadotInput;
	Eigen::Vector4d ThetaddotInput;
	Eigen::Vector2d FeedFwd;
	jt_type dynFeedFWD;

	virtual void operate() {
		tmp_theta_pos = this->jpInputDynamics.getValue();
		ThetaInput << tmp_theta_pos[0], tmp_theta_pos[1], tmp_theta_pos[2], tmp_theta_pos[3];
		Y = calculate_Y_matrix(ThetaInput);
		P = initialize_pi();
		
		FeedFwd = Y * P;
	
		dynFeedFWD[1] = FeedFwd[0];
		dynFeedFWD[3] = FeedFwd[1];

		this->dynamicsFeedFWDValue->setData(&dynFeedFWD);

	}
private:
	DISALLOW_COPY_AND_ASSIGN(Dynamics);
};
#endif /* DYNAMICS_H_ */