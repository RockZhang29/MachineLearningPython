from socket import *
from time import ctime
from datetime import datetime
from typing import Sequence
from egm_pb2 import EgmRobot, EgmHeader, EgmFeedBack, EgmPlanned, EgmMotorState, EgmMCIState, EgmTestSignals, EgmRapidCtrlExecState, EgmMeasuredForce, \
                    EgmSensor, EgmSpeedRef, EgmJoints, EgmCartesianSpeed, \
                    EgmSensorPathCorr, EgmPathCorr, EgmCartesian
from google.protobuf.json_format import MessageToDict

HOST = '127.0.0.1'
PORT = 6510
BUFSIZ = 1024
ADDR = (HOST,PORT)

sequence_num = 0

udpSerSock = socket(AF_INET,SOCK_DGRAM)
udpSerSock.bind(ADDR)

while True:
    data, addr = udpSerSock.recvfrom(BUFSIZ)
    # print('Received from %s:%s.' % addr )
    # print(data)
    # egm_robot = EgmRobot()
    # egm_robot.ParseFromString(data)
    # dict_obj = MessageToDict(egm_robot)
    # print(dict_obj)

    egm_sensor_path_corr = EgmSensorPathCorr()
    egm_sensor_path_corr.header.seqno = sequence_num
    egm_sensor_path_corr.header.tm = datetime.now().microsecond
    egm_sensor_path_corr.header.mtype = EgmHeader.MessageType.MSGTYPE_PATH_CORRECTION
    egm_sensor_path_corr.pathCorr.pos.x = 0.0
    egm_sensor_path_corr.pathCorr.pos.y = 10.0
    egm_sensor_path_corr.pathCorr.pos.z = 10.0
    egm_sensor_path_corr.pathCorr.age = 0

    udpSerSock.sendto(egm_sensor_path_corr.SerializeToString(), addr)

    sequence_num += 1



