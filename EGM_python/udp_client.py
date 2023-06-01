import socket
import time

 
host  = '127.0.0.1' # 这是客户端的电脑的ip
port = 12345 #接口选择大于10000的，避免冲突
bufsize = 4096  #定义缓冲大小
 
addr = (host,port) # 元祖形式
udpClient = socket.socket(socket.AF_INET,socket.SOCK_DGRAM) #创建客户端
while True:
    data = input('>>> ')
    if not data:
        break
    data = data.encode(encoding="utf-8") 
    udpClient.sendto(data,addr) # 发送数据
    data,addr = udpClient.recvfrom(bufsize) #接收数据和返回地址
    print(data.decode(encoding="utf-8"),'from',addr)
 
udpClient.close()