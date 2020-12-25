
      SUBROUTINE NEOPT(M,II,INDEX2,GB,AN,
     *RK,HA,RKI,HAI,STAGE,GKA,GDK,GDA,FF)
      DIMENSION GB(10),AN(10,3),RK(20),HA(20),
     *RKI(10),STAGE(24),GKA(4,30),GDK(6,30),GDA(6,30)
      DIMENSION FF(24),HAI(10)
      REAL LA2,LA4,KAG
      COMMON/N/STUP(160),IZY
      COMMON/A1/D1/A2/R1/A3/YY,H1,H2,T1,T2,CK
      COMMON/A5/HRK
      COMMON/A6/R2,DK1,DK2,DK4,DDK1,DDK2,DDK4
      COMMON/A7/UK1/A8/ALFA1/A9/RRR,AT1,UK2,UK4/A10/AKH
      COMMON/A11/C2AO,C4AO,ALFA4O/A12/AKP/A13/AK
      COMMON/A14/HB3,HB4,CA,BT3,BT4
      COMMON/A15/AKR/A43/FGK,FGKA
      COMMON/A16/SIGMAK,SIGMAA,DSITAK,DSITAA
      COMMON/A17/QLA4,ALFA4, QLA1
      COMMON/A18/AA/BLOK/ARP,A1,R,UKP
      COMMON/A21/KRAN(40)
      COMMON/A22/C(4)
      COMMON/A23/INDEX6,INDEX7,INDEX8,INDEX9,GUIDE
      COMMON/A24/IL,IM,IN,KM,KMI
      COMMON/A25/SL1(7),SL2(6),DL1(11),DL2(10)
      COMMON/A27/SIGMA,DSIGMA,DELTQ
      COMMON/A28/BET2,PPI
      COMMON/A30/C1A,AML1
      COMMON/A31/F1,XFK,XFA/A32/ANT,AKG
      COMMON/A33/QLA1C,HQCK,HQC
      COMMON/A34/PM,P1INIT
      COMMON/A35/IJ
      COMMON/A36/IPEP,KKM,SGM(20),DHQMIN(20),INDEXN
      COMMON/A37/WAG1,WAG2,IND8,IND8A
      COMMON/A38/IND8R
      COMMON/A40/LA4
      COMMON/BIS/BIS(30)
      COMMON/C1/III,INDEXC,HQKP,HQKPA
      COMMON/C3/HQK,HQA
      COMMON/STUP/IQPN
      DOUBLE PRECISION STUP,COD
  112 FORMAT('N STUP',4X,4HQLA1,8X,2HPI,7X,3HKPD,8X,2HHT,8X,2HT2,
     17X,3HC1A,8X,2HHQ,6X,4HHQCK,7X,3HHQC,5X,5HHQKPK,4X,5HHQKPA)
  113 FORMAT(2X,I2,1X,11F10.4)
  933 FORMAT(6F20.6)
      IIOUT=0
      PA=57.296
      HQ=0.
      HQCK=0.
      HT=0.
      AKPX=0.
      PI=0.
      C2A=0.
      PIK=0.
      ALFA2=0.
      IND8=0
      IJ=0
      DK=1E4
      DA=2E4
      HQKP=0.
      HQKPA=0.
      HQC=0.
      KAG=GDK(2,II)
      HAG=GDA(2,II)
      DK1=RK(1)               ! 转子进口处外缘直径Dk1
      DK2=RK(2)               ! 转子出口处外缘直径Dk2
      DDK1=RK(3)              ! 转子进口处轮毂比d1
      DDK2=RK(4)              ! 转子出口处轮毂比d2
      DK4=HA(2)               ! 静子出口处机匣直径Dk2
      DDK4=HA(4)              ! 静子出口处轮毂比d2
      D1=RK(3)
      D2=RK(4)
      R1=SQRT((1+D1**2)/2)    ! 入口处工作轮中间相对半径，以面积为基准
      R2=SQRT((1+D2**2)/2)    ! 出口处工作轮中间相对半径，以面积为基准
      HRK=RK(9)               ! 转子叶片展弦比
      XFK=RK(11)              ! 转子叶片中线的最大挠度相对坐标
      XFA=HA(11)              ! 静子叶片中线的最大挠度相对坐标
      CK=RK(12)               ! 转子叶片型面最大相对厚度
      CA=HA(12)               ! 静子叶片型面最大相对厚度
      YY=RK(8)                ! 工作轮叶片弦长
      H1=RK(9)                ! 工作轮展弦比
      HB3=HA(9)               ! 导向器入口展弦比
      GG1=RK(7)               ! 工作轮叶栅稠度
      GG3=GG1*R1*DK1/(R2*DK2) ! 导向器入口叶栅稠度
      H2=DK2*(1-D2)/2         ! 转子出口流通半径
      T1=YY/GG1               ! 工作轮入口节距
      T2=YY/GG3               ! 导向器入口节距
      HB4=HB3*DK4*(1-DDK4)/(DK2*(1-DDK2))   ! 静子出口展弦比
      BT3=HA(7)               ! 静叶入口叶栅稠度
      BT4=BT3*SQRT((1.+DDK2**2.)/
     *(1.+DDK4**2.))*DK2/DK4  ! 静子出口叶栅稠度
      IF(RK(16).EQ.0) THEN
          UPR2=0.1
      ELSE
          UPR2=2
      END IF
      AT1=GB(2)                           ! 入口气体总温
      RRR=GB(6)                           ! 气体常数R
      ALFA1=GB(4)                         ! 工作轮进口处气流角
      ALFA1=ALFA1/PA                      ! 换算成弧度制
      A1=ALFA1                            ! 入口气流角
      F1=3.14159*DK1**2.*(1.-D1**2.)/4.   ! 转子进口流通面积
      AM=SQRT(AK*(2/(AK+1))**((AK+1)/(AK-1)))      !求流量系数K=AM/SQRT(R)
      IF(II.LE.IL) AA=AM*F1/SQRT(AT1)*SIN(ALFA1)
      IF(II.GT.1) AA=AA*(1.-BIS(II))      ! 级间抽气
      UK2=UK1*DK2/DK1
      UK4=UK1*DK4/DK1
      QLA=AN(M,2)         ! 流量系数
      ANT=AN(M,1)         ! 转速百分比
      AKG=GKA(2,II)       ! 空气流量储存系数
      AKH=GKA(3,II)       ! 理论压头减少系数
      IF(KRAN(5).EQ.1) THEN
      WRITE(7,933) DK1,DK2,DK4,D1,D2,DDK1,DDK2,
     *DDK4,R1,R2,HRK,XFK,XFA,CA,CK,YY,T1,T2,
     *H1,HB3,HB4,BT3,BT4,UPR2,AT1,ALFA1,RRR,AK,AKP,AKH,AKG,UK1,UK2,UK4
      WRITE(7,933) QLA,ANT,UK1
      END IF
      ! 从已有数据中读取以下参数
      C1AO=STAGE(1)
      HTO=STAGE(5)
      AML1=STAGE(6)
      C2AO=STAGE(9)
      ADST=STAGE(4)
      B2O=STAGE(15)
      C4AO=STAGE(19)
      ALFA4O=STAGE(2)
      RE=STAGE(17)        ! 雷诺数?
      PPI=STAGE(22)       ! 总压比
      ! 计算当前级速度及流量参数
      QLA1=QLA*AKG*AA*SQRT(AT1)/(AM*F1*PPI*SIN(ALFA1))    ! 级入口流量系数
      AKR=SQRT(2*9.81*RRR*AT1*AK/(AK+1))                  ! 级入口处临界声速
      AL1=RLQMDA(QLA1)                                    ! 入口处无因次速度数
      C1A=AL1*AKR*SIN(ALFA1)/UK1                          ! 转子入口相对流量系数
      HQ=C1A/C1AO                                         ! 入口处速度系数与最佳速度之比
      ! 计算临界参数
      CALL CRITIC(II,HQ,QLA1,C1AO,B2O,HTO,ADST,KAG,HAG,RE,
     *UPR2,M,RK,HA,RKI,HAI,STAGE,AN,GB,GKA,
     *INDEX2,INDEX3,C1U,C2A,BETA1,BETA2,ALFA2,LA2,AT2,LA4,
     *C4A,PI,PIK,PIT,HT,AKPX)         ! 临界？
      IF(IJ.EQ.0) GO TO 250
 1111 FORMAT('      环形通道阻塞面积')
      IF(KRAN(39).EQ.1) WRITE(7,1111)
      IM=100
      KM=2
      INDEX2=1
      III=II
      GO TO 210
  250 CONTINUE
      IF(INDEX6.EQ.1) GO TO 260
      IF(III.NE.II.OR.ABS(DELTQ).GT.GB(9)) GO TO 251
  260 IF(III.NE.II.OR.ABS(DSIGMA)*INDEX6.GT.GB(10)) GOTO 251
      IF(IM.EQ.100) INDEXN=III
      III=0
      GO TO 201
  251 IF(INDEX6.EQ.0) GO TO 252
      IF(II.NE.IPEP) GO TO 210
      IM=IPEP
      KM=SGM(IPEP)
      IPEP=0
      GO TO 201
  252 IF(IND8.EQ.0) GO TO 210
      GB(9)=1E-1*GB(9)
      INDEX7=0
      DO J=1,10
          DL2(J)=0
          DL1(J)=0
      END DO
      DL1(11)=0
  210 CONTINUE
      IF(INDEX2.EQ.1) RETURN
      IF(IM.NE.II) GO TO 201
      IF(KM.EQ.1.AND.HQK.GT.0.5) GO TO 201
      IF(KM.EQ.2.AND.HQA.GT.0.5) GO TO 201
      IF(INDEX6.EQ.1) GO TO 204
      IF(GB(9).LE.1E-6) GO TO 202
      IF(INDEX8.NE.1) GO TO 201
      DELTQ=0.25*GB(9)
      GB(9)=0.1*GB(9)
      GUIDE=0.
      INDEX7=0
      INDEX8=0
      IM=0
      KM=0
      GO TO 201
  202 INDEXN=II
      GO TO 201
  204 IF(ABS(DSIGMA).GT.GB(10)) GO TO 201
      IF(GB(10).LE.1E-6) GO TO 202
      DSIGMA=0.25*GB(10)      ! 提高精度
      GB(10)=0.1*GB(10)       ! 提高精度
      INDEX9=0
      IM=0
      KM=0
  201 CONTINUE                ! 满足?要求
      IF(HQKP.GT.0.5) DK=HQKP
      DK=DK-HQ
      IF(HQKPA.GT.0.5) DA=HQKPA
      DA=DA-HQ
      DHQMIN(II)=DK
      SGM(II)=SIGMAK+1.
      IF(DHQMIN(II).GT.DA) DHQMIN(II)=DA
      IF(DHQMIN(II).EQ.DA) SGM(II)=SIGMAA+2.
      CALL BRIDGE(II,STAGE,FF,GB)     ! ??????
      IF(HQ.LE.HQCK) INDEXC=II
      GB(1)=PI*GB(1)
      GB(2)=AT2
      GB(4)=ALFA4
      PPI=PPI*PI
      STAGE(22)=PPI
      STUP(IZY)=COD(HQ,4,STAGE(1),4,5)
      STUP(IZY+1)=COD(HQC,4,AKPX,4,6)
      HZ=HT*GKA(3,II)
      HZO=STAGE(5)*GKA(3,II)
      STUP(IZY+2)=COD(HZ,4,HZO,4,5)
      STUP(IZY+4)=COD(HQCK,4,ALFA1*PA,3,6)
      IF(IQPN.EQ.5) THEN
        STUP(IZY+5)=COD(HQKP,4,HQKPA,4,6)
      ELSE
        STUP(IZY+5)=COD(STAGE(15)*PA,2,STAGE(2)*PA,2,6)
      ENDIF
      STUP(IZY+3)=COD(PI,4,STAGE(4),4,5)
      STUP(IZY+6)=COD(C2A,4,LA2,4,5)
      STUP(IZY+7)=COD(PIK,4,ALFA2*PA,2,5)
      ! ***************** 输出该级参数 **************************************************
  700 FORMAT(10X,I2,F12.6,F16.6,2F10.6,2F12.6)
  701 FORMAT(10X,2HII,6X,2HT1,12X,2HP1,10X,5HALFA1,
     *6X,4HQLA1,8X,1HQ)
      IF(II.EQ.1.AND.IIOUT.EQ.1.) WRITE(7,701)
      AMK=0.0404
      QQQ=QLA1*F1/1000**2*AMK*P1INIT*PPI/PI/SQRT(AT1)/AKG*SIN(ALFA1)
      IF(IIOUT.EQ.1.) WRITE(7,700)II,AT1,P1INIT*PPI,ALFA1,QLA1,QQQ
      ! ********************************************************************************
	RETURN
      BETA2=BETA(C2A,R2,ALFA2)
      WRITE(7,113)II,QLA1,PI,AKPX,(beta1-rk(5))*180/3.14,
     *(beta1-rk(5))*180/3.14,(ha(6)-alfa4)*180/3.14
      RETURN
      END