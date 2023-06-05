from socket import *
from abb_robot_client.egm import EGM

egm = EGM()
while True:
    res, RS = egm.receive_from_robot(timeout=0.1)
    if res and RS.rapid_running:
        pos = RS.cartesian[0]
        tm = RS.robot_message.header.tm
        print(f"egm 读取的位置：{pos}，时间：{tm}")

