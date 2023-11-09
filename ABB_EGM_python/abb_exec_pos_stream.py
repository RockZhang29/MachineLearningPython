from abb_robot_client.egm import EGM
import abb_motion_program_exec as abb
import time
import numpy as np

sensor_frame = abb.pose([0,0,0],[1,0,0,0])
egm_config = abb.EGMPathCorrectionConfig(sensor_frame)
mp = abb.MotionProgram(egm_config = egm_config)

egm = EGM()
x,y,z = [],[],[]
count=0
t1 = time.perf_counter()
t2=t1
while (t2 - t1) < 20: 
    t2 = time.perf_counter() 
    res, RobotState = egm.receive_from_robot(timeout=0.1)
    if res:
        if RobotState.rapid_running:
            pos = RobotState.cartesian[0] # ndarray,[x,y,z]
            tm=RobotState.robot_message.header.tm  # 读取时间戳
            count=count+1
            print(f"egm 读取的位置：{pos},时间：{tm}")
    
print(f"egm stram读取的点位数量: {count}")

        
        

        
        
