from socket import *
from time import ctime
from datetime import datetime
from abb_robot_client._egm_protobuf.egm_pb2 import EgmRobot, EgmHeader, EgmFeedBack, EgmPlanned, EgmMotorState, EgmMCIState, EgmTestSignals, EgmRapidCtrlExecState, \
                    EgmSensor, EgmSpeedRef, EgmJoints, EgmCartesianSpeed, \
                    EgmSensorPathCorr, EgmPathCorr, EgmCartesian
from google.protobuf.json_format import MessageToDict

HOST = '127.0.0.1'
PORT = 6510
BUFSIZ = 1024
ADDR = (HOST,PORT)

udpSerSock = socket(AF_INET,SOCK_DGRAM)
udpSerSock.bind(ADDR)

while True:
    data, addr = udpSerSock.recvfrom(BUFSIZ)
    egm_robot = EgmRobot()
    egm_robot.ParseFromString(data)
    dict_obj = MessageToDict(egm_robot)
    # print(dict_obj)
    if egm_robot.planned.cartesian.pos.x<700:
        egm_sensor = EgmSensor()
        egm_sensor.header.seqno = 0
        egm_sensor.header.tm = datetime.now().microsecond
        egm_sensor.header.mtype = EgmHeader.MessageType.MSGTYPE_CORRECTION
        egm_sensor.planned.cartesian.pos.x = egm_robot.planned.cartesian.pos.x + 10
        egm_sensor.planned.cartesian.pos.y = egm_robot.planned.cartesian.pos.y
        egm_sensor.planned.cartesian.pos.z = egm_robot.planned.cartesian.pos.z
        egm_sensor.planned.cartesian.orient.u0 = egm_robot.planned.cartesian.orient.u0
        egm_sensor.planned.cartesian.orient.u1 = egm_robot.planned.cartesian.orient.u1
        egm_sensor.planned.cartesian.orient.u2 = egm_robot.planned.cartesian.orient.u2
        egm_sensor.planned.cartesian.orient.u3 = egm_robot.planned.cartesian.orient.u3

        udpSerSock.sendto(egm_sensor.SerializeToString(), addr)
    else:
        print("超出限位，立刻停止")
        break


