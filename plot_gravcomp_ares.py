import matplotlib.pyplot as plt
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("file_name")
args = parser.parse_args()


# Assuming the data is stored in a text file named 'output.txt'
file_path = args.file_name

# Lists to hold the parsed data for only the 2nd and 4th joints
time_log = []
jp_current_2 = []
jp_current_4 = []
jt_sum_2 = []
jt_sum_4 = []
gravity_2 = []
gravity_4 = []
dynamics_feed_fwd_2 = []
dynamics_feed_fwd_4 = []

# Reading and parsing the data from the file
with open(file_path, 'r') as file:
    for line in file:
        # Assuming the tuple is written in a format like: 
        # (time, jp_current, jp_desired, jt_sum, gravity, control_jt, dynamics_feed_fwd)
        data = eval(line.strip())
        time_log.append(data[0])
        # Append the data for the 2nd and 4th joints
        jp_current_2.append(data[2])  # 2nd joint (index 1)
        jp_current_4.append(data[4])  # 4th joint (index 3)
        jt_sum_2.append(data[6])  # 2nd joint (index 1)
        jt_sum_4.append(data[8])  # 4th joint (index 3)
        gravity_2.append(data[10])
        gravity_4.append(data[12])
        dynamics_feed_fwd_2.append(data[14])
        dynamics_feed_fwd_4.append(data[16])

# Convert lists to numpy arrays for easier handling
time_log = np.array(time_log)
jp_current_2 = np.array(jp_current_2)
jp_current_4 = np.array(jp_current_4)
jt_sum_2 = np.array(jt_sum_2)
jt_sum_4 = np.array(jt_sum_4)
gravity_2 = np.array(gravity_2)
gravity_4 = np.array(gravity_4)
dynamics_feed_fwd_2 = np.array(dynamics_feed_fwd_2)
dynamics_feed_fwd_4 = np.array(dynamics_feed_fwd_4)

# Plotting the data
# Plotting the data
plt.figure(figsize=(12, 16))

# 1. Position for the 2nd joint
plt.subplot(4, 2, 1)
plt.plot(time_log, jp_current_2, label='Position 2nd Joint')
plt.xlabel('Time')
plt.ylabel('Position [rad]')
plt.title('2nd Joint Position')
plt.legend()

# 2. Position for the 4th joint
plt.subplot(4, 2, 2)
plt.plot(time_log, jp_current_4, label='Position 4th Joint')
plt.xlabel('Time')
plt.ylabel('Position [rad]')
plt.title('4th Joint Position')
plt.legend()

# 3. Torque sum for 2nd joint
plt.subplot(4, 2, 3)
plt.plot(time_log, jt_sum_2, label='Torque Sum 2nd Joint')
plt.xlabel('Time')
plt.ylabel('Torque [Nm]')
plt.title('2nd Joint Torque Sum')
plt.legend()

# 4. Torque sum for 4th joint
plt.subplot(4, 2, 4)
plt.plot(time_log, jt_sum_4, label='Torque Sum 4th Joint')
plt.xlabel('Time')
plt.ylabel('Torque [Nm]')
plt.title('4th Joint Torque Sum')
plt.legend()

# 5. Gravity vs Dynamics Feed Fwd for 2nd joint
plt.subplot(4, 2, 5)
plt.plot(time_log, gravity_2, label='Gravity 2nd Joint')
plt.plot(time_log, dynamics_feed_fwd_2, label='Dynamics Fwd 2nd Joint')
plt.xlabel('Time')
plt.ylabel('Torque [Nm]')
plt.title('2nd Joint: Gravity vs Dynamics Fwd')
plt.legend()

# 6. Gravity vs Dynamics Feed Fwd for 4th joint
plt.subplot(4, 2, 6)
plt.plot(time_log, gravity_4, label='Gravity 4th Joint')
plt.plot(time_log, dynamics_feed_fwd_4, label='Dynamics Fwd 4th Joint')
plt.xlabel('Time')
plt.ylabel('Torque [Nm]')
plt.title('4th Joint: Gravity vs Dynamics Fwd')
plt.legend()

# Optional: compare torque sum and dynamics feedforward (extra plots)
plt.subplot(4, 2, 7)
plt.plot(time_log, jt_sum_2, label='Torque Sum 2nd')
plt.plot(time_log, dynamics_feed_fwd_2, label='Dyn Fwd 2nd')
plt.xlabel('Time')
plt.ylabel('Torque [Nm]')
plt.title('2nd Joint: Torque Sum vs Dynamics Fwd')
plt.legend()

plt.subplot(4, 2, 8)
plt.plot(time_log, jt_sum_4, label='Torque Sum 4th')
plt.plot(time_log, dynamics_feed_fwd_4, label='Dyn Fwd 4th')
plt.xlabel('Time')
plt.ylabel('Torque [Nm]')
plt.title('4th Joint: Torque Sum vs Dynamics Fwd')
plt.legend()

plt.tight_layout()
plt.show()
