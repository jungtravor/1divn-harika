      SUBROUTINE OPTIMA(II,RK,HA,RKI,HAI,
     *GB,GKA,X,Y,GDK,GDA,STAGE,INDEX1,INDEX4,OPTIMS,GKI,M)
      ! 确定级的最佳点
      DIMENSION OPTIMS(10),GKI(10)
      DIMENSION RK(20),HA(20),RKI(10),HAI(10),GB(10),GKA(4,30),
     *X(16),Y(16),GDK(6,30),GDA(6,30),GD(6),GK(4),STAGE(24)
      COMMON/A22/C(4)
      COMMON/A20/ARKI(30,10),AHAI(30,11)
      COMMON/A7/UK1/A8/ALFA1/A9/RRR,AT1,UK2,UK4/A10/AKH
      COMMON/A16/SIGMAK,SIGMAA,DSITAK,DSITAA
      COMMON/BLOK/AP1,A1,R,UKR
      DO JJ=1,10
          HAI(JJ)=0.
          RKI(JJ)=0.
      END DO
      DO K=1,4
          GK(K)=GKA(K,II)
      END DO
      IF(INDEX4.LT.1) THEN
          CALL GEOM (RK,X,Y,GD,HA(19))    ! 计算动叶栅几何参数
          DO K=1,6
              GDK(K,II)=GD(K)
          END DO
          CALL GEOM(HA,X,Y,GD,HA(20))     ! 计算静叶栅几何参数
          DO K=1,6
              GDA(K,II)=GD(K)
          END DO
      END IF
      UK2=UK1*RK(2)/RK(1)             ! 动叶进口速度
      UK4=UK1*HA(2)/RK(1)             ! 静叶进口速度
      CALL POINT(RK,HA,GK,GB,GDK,
     *GDA,HAI,RKI,STAGE,II,INDEX1,OPTIMS,M)
      IF(ABS(GK(1)-GKA(1,II)).LT.0.01) GO TO 70
      IF(GK(1).LT.GKI(II)) GKI(II)=GK(1)
   70 CONTINUE
   71  FORMAT(8F9.4)
      RETURN
      END
      
      SUBROUTINE DELTAK(RK,GK,J,E,C2C1,B20,DELT,BETA2)
      ! 落后角及出口相对气流角的确定
      DIMENSION RK(20),GK(4)
      P=57.296
      B20=B20*P
      RK(6)=RK(6)*P
      E=E*P
      IF(J-1)1,1,2
    2 IF(RK(16))10,11,10
   10 XX=1.
    6 XC=0.455/(C2C1*RK(2)/RK(1)-0.545)
      GOTO 7
   11 XX=1.1
      IF(RK(3)-0.4)3,3,4    !RK(3)相对轮毂直径，即轮毂比
    3 XC=1.
      GOTO 7
    4 IF(RK(3)-0.75)5,5,6
    5 XC=1.3*(RK(3)-0.4)/(C2C1*RK(2)/RK(1)-0.545)-2.86*RK(3)+2.145
    7 IF(B20-60.)13,14,14
   13 XB=1.+0.0015*(60.-B20)+0.0003*(60.-B20)**2.
      GOTO15
   14 XB=1.
   15 XD=XC*XX*XB        !落后角修正系数
      IF(GK(1))17,16,17
   17 XD=XD*GK(1)  !落后角校正系数
      GO TO 16
    1 XD=1.
   16 DELT=(0.92*RK(11)**2.-0.002*RK(6)+0.18)*XD/(SQRT(RK(7))/E-0.002*
     *XD)               !落后角!!!!!!!!!!!!      PDF-P16落后角计算公式
      BETA2=RK(6)-DELT
      DELT=DELT/P
      RK(6)=RK(6)/P
      BETA2=BETA2/P   !叶栅出口气流角
      E=E/P
      RETURN
      END
      
      SUBROUTINE AI0(RK,DI,BETA2,AI0K,BETA1,FER)
      ! 最佳入口攻角计算 PDF - P17,2.入口相对气流角及最佳攻角的确定
      DIMENSION RK(20)
      P=57.296
      BETA2=BETA2*P
      RK(5)=RK(5)*P
      DI=DI*P
       BT=RK(7)       ! 稠度
      IF(BT.GT.0.55) THEN         ! 根据稠度计算叶形弯角
      E=0.5*(0.649+SQRT(3.2*RK(7)-1.373))
      ELSE
      E=1.155*BT
      END IF
      BE=0.001997*BETA2**2.+0.170433*BETA2+2.884
      AI0K=(RK(7)*((RK(5)/30.)**2.+1.5)+10.*RK(3)*(1.-(BETA2-
     *RK(5))/(E*BE))/SQRT(RK(7))-FER+DI)/(1.+10.*RK(3)
     */(SQRT(RK(7))*E*BE)+0.0022*RK(5)*RK(7))+RK(18)
      AI0K=AI0K/P
      RK(5)=RK(5)/P
      BETA1=(RK(5)-AI0K)
      BETA2=BETA2/P
      RETURN
      END
      
      FUNCTION SPID(R1,BETA1,B1)
      ! 级入口叶尖轴向相对速度的确定
      SPID=R1/(CTG(BETA1)+CTG(B1))
      RETURN
      END
      
      SUBROUTINE AL1(C1A0,R1,B10S,AM10,AL10)
      ! 转子入口气流的马赫数及速度系数的确定
      ! ---   C1A0    入口轴向相对速度
      ! ---   R1      入口相对平均半径
      ! ---   B10S    入口相对气流角
      ! ---   AM10    计算结果：入口气流马赫数
      ! ---   AL10    计算结果：入口气流马赫数
      COMMON/A7/UK1/BLOK/AP1,A1,R,UKP/A13/AK
      COMMON/A9/RRR,AT1,UK2,UK4
      COMMON/T1/T1
      REAL*8 AKK
      T1=AT1-(C1A0*UK1/SIN(A1))**2.*(AK-1)/(2.*AK*9.81*R)
      B10S=1.570796-ATAN(R1/C1A0-CTG(A1))
      AKK=AK*9.81*R*T1
      AM10=C1A0*UK1/(SIN(B10S)*SQRT(AKK))
      AL10=SQRT(((AK+1.)*AM10**2./2.)/(1.+(AK-1.)*AM10**2./2.))
      RETURN
      END
      
      FUNCTION AKPD(RK,AL10)
      ! 转子效率计算
      DIMENSION RK(20)
      XQ=AL10
      XR=RK(3)
      IF(RK(16))1,2,1
      ! 亚音速计算
    2 IF(RK(3)-0.55)3,3,4
    3 RK(3)=0.55
      GOTO 5
    4 IF(RK(3)-0.7)5,6,6
    6 RK(3)=0.7
    5 IF(AL10-0.72)7,7,8
    7 IF(AL10.LT.0.62) AL10=0.62
      ! 亚音速叶形，入口最佳速度系数小于0.62
      BB=0.222*RK(3)-0.255
      CC=0.0267*RK(3)+0.9573
      AKPD=BB*AL10+CC
      GO TO 16
    8 AA=1.14*RK(3)-3.752
      BB=5.238-1.47*RK(3)
      CC=0.653*RK(3)-1.053
      AKPD=AA*AL10**2.+BB*AL10+CC
      GO TO 16
      ! 超音速计算
    1 IF(AL10-0.6)9,9,10
    9 XKPD=0.89
      GOTO13
   10 XKPD=0.905-0.025*AL10
      IF(AL10.GT.1.1) XKPD=0.8775-(AL10-1.1)**2
   13 IF(AL10-RK(17))14,14,15
   14 AKPD=XKPD
      GO TO 16
   15 AKPD=XKPD-0.25*(AL10/RK(17)-1.)
   16 AL10=XQ
      RK(3)=XR
      RETURN
      END
      
      FUNCTION HP(AHP,AL10)
      ! 根据流道尺寸对效率的修正值（P20,3-23）
      HR=AL10**2*(25.-AHP)/144.
      IF(HR.LT.0) HR=0.
      HP=HR
      RETURN
      END
      
      SUBROUTINE OKPD(RK,C1A0,C1A0K,AL10,ANW,OK,FK)
      ! 效率修正系数Ok计算
      DIMENSION RK(20)
      AL=AL10        !进口无因次速度数
      ANW=C1A0/C1A0K
      IF(RK(16).EQ.0) THEN
          IF(AL10.LE.0.5) AL10=0.5
          FK=AL10**2.-0.97*AL10+0.3
      ELSE
          IF(AL10.LE.0.9) AL10=0.9
          FK=0.46*AL10-0.35
      END IF
      OK=1-10.*(ABS(ANW-1.)**1.334)*FK  !转子效率修正系数
      AL10=AL
      RETURN
      END
      
      SUBROUTINE PPZPKC(HT,AKPX,QLA1,C1A,C2AOA,
     *PIK,ALFA2,LA2,C2A,C2U,AT2,
     *PI,PIT,ALFA4,LA4,C4A,II)
      ! 转子出口参数计算
      REAL LA2,LA4,LA2A
      COMMON/A2/R1/A6/R2,DK1,DK2,DK4,DDK1,DDK2,DDK4
      COMMON/A7/UK1/A8/ALFA1
      COMMON/A9/RRR,AT1,UK2,UK4/A10/AKH
      COMMON/A11/C2AO,C4AO,ALFA4O/A13/AK
      HZ=HT*UK1**2*AKH/9.81
      HF=HT*AKH
      C1U=C1A*COS(ALFA1)/SIN(ALFA1)  !进口周向相对速度系数
      DF2=DK2/DK1
      C2U=(HF+R1*C1U)/(R2*DF2**2)    !出口周向相对速度系数
      ! 反动度
      TAW=1-(C1U+C2U*DF2)/(2.0*R1)
      IF(TAW.LT.0.3) TAW=0.3
      IF(TAW.GT.0.8) TAW=0.8
      ! 转子效率
      AKPDK=1-(1-AKPX)*(0.8*TAW+0.1)
      HAD=HZ*AKPX     !级的等熵功
      HADK=HZ*AKPDK   !转子的等熵功
      FPI=AK*AT1*RRR/(AK-1)
      FPP=AK/(AK-1)
      PI=(1+HAD/FPI)**FPP   !级压比
      PIK=(1+HADK/FPI)**FPP  !转子压比
      PIT=(1+HZ/FPI)**FPP    !级理论压比
      AT2=AT1+HZ/(FPP*RRR)   !出口总温
      SMUL=1.
      IF(HZ.LT.0.) SMUL=-1.
      IF(ABS(PI-1.).LT.1.E-5) THEN
          PI=1.+SMUL*1.E-5
          PIK=1+SMUL*(HADK/HAD)*1.E-5
          PIT=1.+SMUL*(HZ/HAD)*1.E-5
      END IF
      IF(ABS(AT2-AT1).LT.1.E-5*AT1) AT2=AT1*(1.+SMUL*1.E-5)
      AT2A=AT2-(AK-1)*(C2U*UK2)**2/(2*AK*9.81*RRR)
      WW1=DK1**2/DK2**2
      WW3=(AT2/AT2A)**((AK+1)/(2*AK-2))
      WW2=(1-DDK1**2)/(1-DDK2**2)
      WW4=SQRT(AT2/AT1)
      QLA2A=QLA1*SIN(ALFA1)*WW1*WW2*WW3*WW4/PIK  !求q
      LA2A=RLQMDA(QLA2A)   !转子出口轴向无因次速度数
      WW5=SQRT(2*9.81*RRR*AT2A*AK/(AK+1))
      C2A=LA2A*WW5/UK2  !转子出口轴向相对速度系数
      AKR2=SQRT(2*AK*AT2*9.81*RRR/(AK+1))  !转子出口临界音速
      LA2=UK2*SQRT(C2A**2+C2U**2)/AKR2  !转子出口无因次速度数
      TGB1=C2A/(R2-C1U)
      TGA2=C2A/C2U  !出口气流角正切值
      ALFA2=1.570796-ARCCTG(TGA2)  !出口气流角
      ! 计算出口流量及流速
      CALL AEXIT(C2A,AT2,PI,QLA1,
     *C2AOA,ALFA4,LA4,C4A,II)
      RETURN
      END
      
      FUNCTION ATAKA(RK,AI0K,R1,C1A0K,B1,AL10)
      ! 攻角修正系数计算
      DIMENSION RK(20)
      IF(AL10-0.4)1,1,2
    1 X=1.09
      GOTO5
    2 IF(AL10-0.801)3,3,4
    3 X=1.09-0.375*(AL10-0.4)**2.-0.075*(AL10-0.4)
      GOTO5
    4 X=2.53*AL10-1.93*AL10**2.+0.211
    5 ATAKA=RK(5)-AI0K-ARCCTG(R1/(C1A0K*X)-CTG(B1))
      RETURN
      END
      
      FUNCTION ONWK(RK,ANWK,ANWA,TAW,AKPDK,AKPDA,FK,FA)
      ! 轴向相对速度修正系数计算
      DIMENSION RK(20)
      X=FA/FK
      IF(RK(16))1,3,1
    1 X=0.7*X
    3 AO=AKPDA*X*(1.-TAW)/(AKPDK*TAW)
      ONWK=(1.+AO*ANWA/ANWK)/(1.+AO*(ANWA/ANWK)**2.)
      RETURN
      END