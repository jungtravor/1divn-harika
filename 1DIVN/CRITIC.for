﻿
      SUBROUTINE CRITIC(II,HQ,QLA1,C1AO,B2O,HTO,ADST,KAG,HAG,RE,UPR2,M,
     *RK,HA,RKI,HAI,STAGE,AN,GB,GKA,
     *INDEX2,INDEX3,C1U,C2A,BETA1,BETA2,ALFA2,LA2,AT2,
     *LA4,C4A,PI,PIK,PIT,HT,AKPX)
      DIMENSION RK(20),HA(20),RKI(10),HAI(10),STAGE(24)
      DIMENSION GB(10),AN(10,3),GKA(4,30)
      COMMON/PLAN/I,IZ,N,NTB,IKBD
      COMMON/A1/D1/A2/R1/A3/YY,H1,H2,T1,T2,CK
      COMMON/A5/HRK
      COMMON/A6/R2,DK1,DK2,DK4,DDK1,DDK2,DDK4
      COMMON/A7/UK1/A8/ALFA1/A9/RRR,AT1,UK2,UK4/A10/AKH
      COMMON/A11/C2AO,C4AO,ALFA4O/A12/AKP/A13/AK
      COMMON/A14/HB3,HB4,CA,BT3,BT4
      COMMON/A15/AKR
      COMMON/A16/SIGMAK,SIGMAA,DSITAK,DSITAA
      COMMON/A17/QLA4,ALFA4,QQLA1
      COMMON/A21/KRAN(40)
      COMMON/A22/C(4)
      COMMON/A23/INDEX6,INDEX7,INDEX8,INDEX9,GUIDE
      COMMON/A24/IL,IM,IN,KM,KMI
      COMMON/A25/SL1(7),SL2(6),DL1(11),DL2(10)
      COMMON/A25V/SIGF
      COMMON/A27/SIGMA,DSIGMA,DELTQ
      COMMON/A28/BET2,PPI
      COMMON/A30/C1A,AML1
      COMMON/A35/IJ
      COMMON/C1/III,INDEXC,HQKP,HQKPA
      COMMON/C2/KSI,SL21(6)/C2R/SL22(6)
      COMMON/C3/HQK,HQA
      COMMON/A31/F1,XFK,XFA/A32/ANT,AKG
       REAL LA2,LA4,KAG
      DIMENSION SL11(6)
      DIMENSION PSL1(7)
  413 FORMAT('   II IL IM IN KM KMI
     *INDEX2 INDEX3 INDEX6 INDEX7 INDEX8 INDEX9')
  414 FORMAT('   GUIDE HQ KAG HAG F1 AT1 ALFA1
     *SIGMA DSIGMA DELTQ')
  700 FORMAT(6F20.6)
  800 FORMAT('   CRITIC   崁梹')
  801 FORMAT (5X,7H HQKPK=,F9.6)
  802 FORMAT (5X,7H HQKPA=,F9.6)
  900 FORMAT('   CRITIC   妿崡厤')
 1000 FORMAT(6I20)
 1975 FORMAT(12X,2X,2HI=,I2,2X,7HICOUNT=,I2,2X,I3,2X,5HQLA1=,F8.6,2X,2HQ
     *=F16.6)
      IOUT=0
      IOFILE=0
      IF(KRAN(35).EQ.1) WRITE(6,800)
      ICOUNT=0
      QLAM=AN(M,2)
      Z=0.
      HQK=0.
      HQA=0.
      IR=0
      DS=DSIGMA
      S=SIGMA
      HTM=3.
      PSL1(1)=0.
      C2AOA=HAI(9)
      HQP=HQ          ! 入口轴向速度与最佳入口轴向速度之比
      IND4=0
      IND5=0
      INDEX3=0
   28 CONTINUE
      ICOUNT=ICOUNT+1
      C1A=HQ*C1AO
      IND1=0
      IND2=0
      IND3=0
      K=0
      KK=0
      QLA1=QQL(C1A)       ! 入口流量系数
      IF(IOUT.EQ.1) WRITE(IOFILE,1975)II,ICOUNT,28,QLA1,CQ(QLA1)
      IF(KRAN(35).EQ.1) WRITE(6,413)
      IF(KRAN(35).EQ.1) WRITE(6,*) '  ICOUNT= ',ICOUNT,' Z= ',Z
      IF(KRAN(35).EQ.1) WRITE(6,1000) II,IL,IM,IN,KM,KMI,
     *INDEX2,INDEX3,INDEX6,INDEX7,INDEX8,INDEX9,
     *IND1,IND2,IND3,IND4,IND5
      IF(KRAN(35).EQ.1) WRITE(6,414)
      IF(KRAN(35).EQ.1) WRITE(6,700) GUIDE,HQ,KAG,HAG,F1,
     *AT1,ALFA1,SIGMA,DSIGMA,DELTQ,C1A,STAGE(23)
      IF(KRAN(35).EQ.1) WRITE(6,700) Z
      IF(IN.EQ.II) GO TO 100
      IF(INDEX3.NE.5) GO TO 201
      IND1=1
      GO TO 49
  201 IF(INDEX3.NE.3) GO TO 30
      HQKP=SL11(1)
      C1AKP=SL11(2)
      AL11KP=SL11(3)
      BET2=SL11(4)
      PITKP=SL11(5)
      PIKKP=SL11(6)
      IND1=1
      GO TO 30
  100 IF(KMI.EQ.1) GO TO 51
      GO TO 101
   30 C1A=C1A/STAGE(23)
      BETA1=BETA(C1A,R1,ALFA1)   !!!!!!! 进口气流角
      TC1=AT1-(C1A*UK1/SIN(ALFA1))**2*(AK-1.)/(2.*9.81*AK*RRR)    !进口静温
      AM11=C1A*UK1/(SIN(BETA1)*SQRT(9.81*AK*RRR*TC1))             !动叶进口相对马赫数
      AL11=AM11*SQRT((AK+1.)/(2.+(AK-1.)*AM11*AM11))              !无因次速度数
      C1A=C1A*STAGE(23)
      IF(INDEX3.EQ.3) GO TO 301
      CALL DOWN1(II,RK,C,C1A,KAG,BETA1,AM11,AMMAX,KM,IM,KK)
      IF(KRAN(35).EQ.1) WRITE(6,1000) KM,IM,KK
      GO TO (68,69),KK     !判断是否大于最大马赫数，垂直段起点
   68 CONTINUE
      IF(RK(16).EQ.0) GO TO 29
      IF((C1A-C(2)).LE.0) GO TO 32
      GO TO 33
   29 CALL ELEVEL(KAG,BETA1,AMKPK)
   31 IF(AM11-AMKPK)32,32,33    !是否大于临界马赫数
   32 IF(INDEX3.EQ.0) GO TO 34
      IF(INDEX3.NE.2) GO TO 40
   34 CALL BLOCK(II,HQ,QLA1,C1AO,B2O,HTO,ADST,RE,UPR2,
     *RK,HA,RKI,HAI,STAGE,
     *ALFA2,C2A,AT2,PI,PIK,PIT,LA4,C4A,HT,AKPX)  ! 各区域参数计算？分离区参数计算
      IF(IOUT.EQ.1) WRITE(IOFILE,1975)II,ICOUNT,34,QLA1,CQ(QLA1)
      LA2=C2A*UK2/(SIN(ALFA2)*SQRT(2.*9.81*RRR*AK*AT2/(AK+1.)))
      CALL DOWN2(II,GB,QLAM,LA2,ALFA2,HAG,AK,AML2,AMMAX,KM,IM,KK)
      GO TO (70,69),KK
   70 CONTINUE
      CALL ALEVEL(HAG,ALFA2,LA2,AK,AMKPA,AM2,K,Z)
      IF(KRAN(35).EQ.1) WRITE(6,700) Z
      IF(K.GT.2.AND.HTM.GE.2.) HTM=HT
      GO TO (35,35,36),K
   35 IF(INDEX3.NE.0) GOTO 39
      SIGMA=S
      DSIGMA=DS
      RETURN
   36 IF(INDEX3.NE.0) GO TO 37
      SL2(1)=QLA1
      SL2(2)=C2A
      SL2(3)=AT2
      SL2(4)=PIK
      SL2(5)=PIT
      SL2(6)=PI
   37 IND2=2
   39 IF(Z.GT.1E-4.AND.ICOUNT.LT.25) GO TO 38
      GO TO 52
   33 IF(INDEX3.NE.0) GO TO 41
  301 CONTINUE
      SL1(1)=C1A
      SL1(2)=BETA1
      SL1(3)=AL11
      CALL BLOCK(II,HQ,QLA1,C1AO,B2O,HTO,ADST,RE,UPR2,
     *RK,HA,RKI,HAI,STAGE,
     *ALFA2,C2A,AT2,PI,PIK,PIT,LA4,C4A,HT,AKPX)
      IF(IOUT.EQ.1) WRITE(IOFILE,1975)II,ICOUNT,301,QLA1,CQ(QLA1)
      LA2=C2A*UK2/(SIN(ALFA2)*SQRT(2.*9.81*RRR*AK*AT2/(AK+1.)))
      SL1(4)=PI
      SL1(5)=PIK
      SL1(6)=PIT
      SL1(7)=LA2
      IF(INDEX3.NE.0) GO TO 330
      DO J=1,7
          PSL1(J)=SL1(J)
      END DO
  330 CONTINUE
      IF(INDEX3.EQ.3) GO TO 51
   41 IND1=1
   40 IF(RK(16).GT.0) GO TO 66
      Z=ABS(AM11-AMKPK)
   65 CONTINUE
      IF(KRAN(35).EQ.1) WRITE(6,700) Z
      IF(Z.GT.1E-4.AND.ICOUNT.LT.25) GO TO 38
      GO TO 67
   66 Z=ABS(C1A-C(2))
      GO TO 65
   67 HQKP=HQ
      C1AKP=C1A
      AL11KP=AL11
      QLA1=QQL(C1A)
      IF(IOUT.EQ.1) WRITE(IOFILE,1975)II,ICOUNT,671,QLA1,CQ(QLA1)
      CALL BLOCK(II,HQ,QLA1,C1AO,B2O,HTO,ADST,RE,UPR2,
     *RK,HA,RKI,HAI,STAGE,
     *ALFA2,C2A,AT2,PI,PIK,PIT,LA4,C4A,HT,AKPX)
      BETA2=BETA(C2A,R2,ALFA2)
      IF(IOUT.EQ.1) WRITE(IOFILE,1975)II,ICOUNT,672,QLA1,CQ(QLA1)
      BET2=BETA2
      PITKP=PIT
      PIKKP=PIK
      SL11(1)=HQKP
      SL11(2)=C1AKP
      SL11(3)=AL11KP
      SL11(4)=BET2
      SL11(5)=PITKP
      SL11(6)=PIKKP
   51 CALL SLOPE1(SL1,II,GKA,C1AKP,BET2,AL11KP,PITKP,PIKKP,F1,
     *C1U,C2A,C2,BETA2,ALFA2,LA2,AT2,ALFA4,LA4,C4A,PI,PIK,PIT,HT,AKPX)
      IF(IJ.GT.0) RETURN
      HQK=1.
      FCA=SL1 (1)
      QLA1=QQL(FCA)
      IF(IOUT.EQ.1) WRITE(IOFILE,1975)II,ICOUNT,51,QLA1,CQ(QLA1)
      IF(HTM.GE.2) HTM=HT
      CALL DOWN2(II,GB,QLAM,LA2,ALFA2,HAG,AK,AML2,AMMAX,KM,IM,KK)
      GO TO (71,69),KK
   71 IF(INDEX3.NE.1) GO TO 220
      IND1=1
      IND4=0
  220 CONTINUE
      IF(INDEX3.GT.2) GO TO 42
      IF(IR.EQ.1) GOTO 42
      SL2(1)=QLA1
      SL2(2)=C2A
      SL2(3)=AT2
      SL2(4)=PIK
      SL2(5)=PIT
      SL2(6)=PI
   42 CALL ALEVEL(HAG,ALFA2,LA2,AK,AMKPA,AM2,K,Z)
      IF(II.NE.IN) GO TO 250
      IF(IR.EQ.0.AND.K.LE.2) GO TO 251
      IR=1
      IF(K.LE.2.AND.ABS(DSIGMA).LT.GB(10)) GO TO 52
  250 CONTINUE
      IF(IND4.EQ.1.AND.K.GT.2) IND2=2
      GO TO (43,43,44),K
   43 IF(IND5.EQ.1) GO TO 102
      IF(ABS(DHQ).LT.0.005) GO TO 300
      IF(INDEX3.EQ.3) IND2=2
  300 CONTINUE
      IF(INDEX3.NE.1) GO TO 44
  251 CONTINUE
      CALL AEXIT(C2A,AT2,PI,QLA1,C2AOA,ALFA4,LA4,C4A,II)
      IF(IOUT.EQ.1) WRITE(IOFILE,1975)II,ICOUNT,251,QLA1,CQ(QLA1)
      GO TO 210
  102 IND4=1
      DSIGMA=-ABS(DSIGMA)
      GO TO  105
  103 IF(ABS(DSIGMA).LT.GB(10)/2.) DSIGMA=ABS(DSIGMA)*2.
  502 SIGMA=SIGMA+DSIGMA
      GO TO 51
   44 IF(KMI.NE.1) GO TO 104
      IF(IN.NE.II) GO TO 104
      IF(KRAN(35).EQ.1) WRITE(6,700) KSI
      IF(KSI.EQ.0) GO TO 450
      C2AKP=SL21(1)
      PIKKP=SL21(2)
      PITKP=SL21(3)
      PIKP=SL21(4)
      AL2KP=SL21(5)
      GO TO 101
  450 CONTINUE
      IF(IR.EQ.1) GO TO 252
      IF(Z.LE.1E-4.OR.ICOUNT.GE.25) GO TO 101
  252 DSIGMA=ABS(DSIGMA)
      IND5=1
      IF(IND4.EQ.0) GO TO 502
  105 DSIGMA=DSIGMA/2.
      GO TO 103
  104 IF(INDEX3.EQ.3) GO TO 47
      HQ=SL11(1)
      FCA=SL11(2)
      QLA1=QQL(FCA)
      IF(IOUT.EQ.1) WRITE(IOFILE,1975)II,ICOUNT,104,QLA1,CQ(QLA1)
   49 CALL BLOCK(II,HQ,QLA1,C1AO,B2O,HTO,ADST,RE,UPR2,
     *RK,HA,RKI,HAI,STAGE,
     *ALFA2,C2A,AT2,PI,PIK,PIT,LA4,C4A,HT,AKPX)
      IF(IOUT.EQ.1) WRITE(IOFILE,1975)II,ICOUNT,49,QLA1,CQ(QLA1)
      LA2=C2A*UK2/(SIN(ALFA2)*SQRT(2.*9.81*RRR*AK*AT2/(AK+1.)))
      CALL ALEVEL(HAG,ALFA2,LA2,AK,AMKPA,AM2,K,Z)
      IF(KRAN(35).EQ.1) WRITE(6,700) Z
      GO TO (45,52,46),K
   46 IND2=4
   47 IF(Z.GT.1.E-4.AND.ICOUNT.LT.25) GO TO 38
   52 C2AKP=C2A
      PIKKP=PIK
      PITKP=PIT
      PIKP=PI
      AL2KP=LA2
      HQKPA=C1A/C1AO
  101 CALL SLOPE2(II,SL2,C2AKP,PIKKP,PITKP,PIKP,AL2KP,C2AOA,
     *PI,PIK,PIT,LA4,C4A,ALFA4,AKPX,HT)
      HQA=1.
      IF(II.NE.IM) GO TO 500
      SL22(1)=C2AKP
      SL22(2)=PIKKP
      SL22(3)=PITKP
      SL22(4)=PIKP
      SL22(5)=AL2KP
       SL22(6)=HQA
      IF(INDEX8.EQ.1.OR.INDEX8.GT.3) GO TO 501
  500 IF(II.NE.IN) GO TO 451
      IF(KMI.NE.1.OR.KSI.GT.0) GO TO 451
  501 KSI=1
      SL21(1)=C2AKP
      SL21(2)=PIKKP
      SL21(3)=PITKP
      SL21(4)=PIKP
      SL21(5)=AL2KP
  451 CONTINUE
 1516 FORMAT('    CRITIC')
 1517 FORMAT(20I4)
 1518 FORMAT(10F12.6)
      IF(HTM.LE.2) HT =HTM
      AT2=SL2(3)
      IF(IR.EQ.1.AND.ABS(DSIGMA).LT.GB(10)) GOTO 211
      IF(PSL1(1).EQ.0) GO TO 210
      DO J=1,7
          SL1(J)=PSL1(J)
      END DO
  210 CONTINUE
  211 IF(II.NE.IN) GO TO 214
      SIGMA=S
      DSIGMA=DS
  214 IF(IM.NE.II) GO TO 212
      IF(INDEX6.EQ.1.AND.ABS(DSIGMA).LT.GB(10)) GO TO 215
      IF(INDEX8.LT.4.AND.INDEX8.NE.1) GOTO 212
  215 CONTINUE
      IF(KM.EQ.1) GO TO 213
      SIGF=PI/PIK
      GO TO 212
  213 SIGF=PIK/PIT
  212 HQ=HQP
      C1A=HQ*C1AO
      QLA1=QQL(C1A)
      IF(IOUT.EQ.1) WRITE(IOFILE,1975)II,ICOUNT,212,QLA1,CQ(QLA1)
      RETURN
   45 IND2=2
      GO TO 47
   38 CALL WAY3(INDEX3,IND1,IND2,IND3,IND4,HQ,DHQ)
      GO TO 28
   69 KK=0
      INDEX2=1
      K=0
      IF(KRAN(35).EQ.1) WRITE(6,1000) II,IL,IM,IN,KM,KMI,
     *INDEX2,INDEX3,INDEX6,INDEX7,INDEX8,INDEX9
      RETURN
      END