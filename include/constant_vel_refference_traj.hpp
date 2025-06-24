#pragma once

#include <barrett/units.h>
#include <barrett/systems.h>

template<size_t DOF>
class constVelRefTrajectory: public systems::System {
	BARRETT_UNITS_TEMPLATE_TYPEDEFS(DOF);

// IO
public:
	Input<double> timef;
	Output<jp_type> referencePTrack;
	Output<jv_type> referenceVTrack;
	Output<ja_type> referenceATrack;

protected:
	typename System::Output<jp_type>::Value* referencePOpValue;
	typename System::Output<jv_type>::Value* referenceVOpValue;
	typename System::Output<ja_type>::Value* referenceAOpValue;

public:
	explicit constVelRefTrajectory(jp_type start_pose, jv_type v_desired, const std::string& sysName = "constRefTrajectory"):
		System(sysName), timef(this), referencePTrack(this, &referencePOpValue), referenceVTrack(this, &referenceVOpValue),
		referenceATrack(this, &referenceAOpValue), start_pos(start_pose), refPTrack(start_pose), v_des(v_desired), refVTrack(v_desired), refATrack(){}

	virtual ~constVelRefTrajectory() { this->mandatoryCleanUp(); }

protected:
	jp_type start_pos, refPTrack;
	jv_type v_des, refVTrack;
	ja_type refATrack;

	virtual void operate() {
		// refATrack.Zero();
		refPTrack = v_des * this->timef.getValue() + start_pos;
		// refVTrack = v_des;

		referencePOpValue->setData(&refPTrack);
		referenceVOpValue->setData(&refVTrack);
		referenceAOpValue->setData(&refATrack);
	}

private:
	DISALLOW_COPY_AND_ASSIGN(constVelRefTrajectory);
};
