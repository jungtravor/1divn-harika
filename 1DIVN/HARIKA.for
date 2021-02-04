      SUBROUTINE HARIKA(GB,GKA,AN,A,B,X,Y,SIG,ALF,NP,IBD,IQP,IWRITE,AKGK
     *,NNS,NNF,IIS,IIF,QQS,QQF,IPE4,KPAH,FKPAH,DSIG,IZAP,LSR,WCP,CU,EP
     *,PR,PB,BRANCH,OPTIMS,GKI,POINTS,BASE,STORE,GAN,PLANE,NOWNV)
      DIMENSION A(20,30),B(20,30),RKI(10),HAI(10),RK(20),HA(20),
     *STAGE(24),GB(10),GKA(4,30),AN(10,3),X(16),Y(16),
     *KPAH(40),FKPAH(10),BRANCH(10),OPTIMS(10)
      DIMENSION FFX(24),AGP(3)
      character*20 numberchar
      COMMON/AKK/AKPX
      ! AKK
      ! ---   AKPX    等熵效率
      COMMON/G/GDK(6,30),GDA(6,30)
      ! G
      ! ---   GDK(6,30)   动叶栅几何参数
      ! ---   GDA(6,30)   静叶栅几何参数
      ! ------    GD?(1,I)    气流折转角
      ! ------    GD?(2,I)    喉部面积
      ! ------    GD?(3,I)    气流折转角
      ! ------    GD?(4,I)    气流折转角
      ! ------    GD?(5,I)    叶栅进出口面积比
      ! ------    GD?(6,I)    进气角与叶型楔角一半的差
      COMMON/IQ/IQP1/ISTUP/IQPN
      COMMON/PLAN/I,IZ,N,NTB,IKBD
      ! PLAN
      ! ---   I       
      ! ---   IZ      压气机总级数
      ! ---   N       特性线数量
      ! ---   NTB     按折合转差率计算第二级的符号???
      ! ---   IKBD    
      COMMON/A1/D1/A2/R1/A3/YY,H1,H2,T1,T2,CK
      ! A1
      ! ---   D1      转子进口处轮毂比d1
      ! A2
      ! ---   R1      入口处工作轮中间相对半径
      ! A3
      ! ---   YY      工作轮叶片弦长
      ! ---   H1      工作轮叶片展弦比？
      ! ---   H2      工作轮出口流道径向高度的一半
      ! ---   T1      工作轮入口节距
      ! ---   T2      导向器入口节距
      ! ---   CK      工作轮型面最大相对厚度
      COMMON/A5/HRK
      ! A5
      ! ---   HRK     工作轮叶片展弦比
      COMMON/A6/R2,DK1,DK2,DK4,DDK1,DDK2,DDK4
      ! A6
      ! ---   DK1     转子进口处机匣直径Dk1
      ! ---   DK2     转子出口处机匣直径Dk2
      ! ---   DK4     级出口处机匣直径Dk4
      ! ---   DDK1    转子进口处轮毂比d1
      ! ---   DDK2    转子出口处轮毂比d2
      ! ---   DDK4    级出口处轮毂比d4
      COMMON/A7/UK1/A8/ALFA1/A9/RRR,AT1,UK2,UK4/A10/AKH
      ! A7
      ! ---   UK1     工作轮进口周向速度
      ! A8
      ! ---   ALFA1   工作轮进口处气流角
      ! A9
      ! ---   RRR     气体常数R
      ! ---   AT1     转子进口总温
      ! ---   UK2     转子出口周向速度
      ! ---   UK4     静子出口周向速度
      ! A10
      ! ---   AKH     理论压头减少系数
      COMMON/A11/C2AO,C4AO,ALFA4O/A12/AKP/A13/AK
      ! A11
      ! ---   C2AO    工作轮出口最佳轴向速度
      ! ---   C4AO    导向器出口最佳轴向速度
      ! ---   ALFA4O  导向器出口处最佳气流角
      ! A12
      ! ---   AKP     压头特性校正系数
      ! A13
      ! ---   AK      等熵过程指数
      COMMON/A14/HB3,HB4,CA,BT3,BT4
      ! A14
      ! ---   HB3     导向器入口处叶片展弦比
      ! ---   HB4     导向器出口处叶片展弦比
      ! ---   CA      导向器叶片型面最大相对厚度
      ! ---   BT3     导向器入口处叶栅稠度
      ! ---   BT4     导向器出口处叶栅稠度
      COMMON/A15/AKR
      ! A15
      ! ---   AKR     级入口处临界声速
      COMMON/A16/SIGMAK,SIGMAA,DSITAK,DSITAA
      ! A16
      ! ---   SIGMAK  ???
      ! ---   SIGMAA  ???
      ! ---   DSITAK  ???
      ! ---   DSITAA  ???
      COMMON/A17/QLA4,ALFA4, QLA1
      ! A17
      ! ---   QLA4    级出口折合流量
      ! ---   ALFA4   静子出口气流角
      ! ---   QLA1    级入口折合流量
      COMMON/A18/AA
      ! A18
      ! ---   AA      ???
      COMMON/A20/ARKI(30,10),AHAI(30,11)
      ! A20
      ! ---   ARKI    ???
      ! ---   AHAI    ???
      COMMON/A21/KRAN(40)
      ! A21
      ! ---   KRAN    ???
      COMMON/A22/C(4)
      ! A21
      ! ---   C       ???
      COMMON/A23/INDEX6,INDEX7,INDEX8,INDEX9,GUIDE
      ! A23
      ! ---   INDEX6  ???
      ! ---   INDEX7  ???
      ! ---   INDEX8  ???
      ! ---   INDEX9  ???
      ! ---   GUIDE   ???
      COMMON/A24/IL,IM,IN,KM,KMI
      ! A24
      ! ---   IL      ???
      ! ---   IM      ???
      ! ---   IN      ???
      ! ---   KM      ???
      ! ---   KMI     ???
      COMMON/A25/SL1(7),SL2(6),DL1(11),DL2(10)
      ! A25
      ! ---   SL1     ???
      ! ---   SL2     ???
      ! ---   DL1     ???
      ! ---   DL2     ???
      COMMON/A25R/PDL1(11),PDL2(10)
      ! A25R
      ! ---   PDL1    ???
      ! ---   PDL2    ???
      COMMON/A25V/SIGF
      ! A25V
      ! ---   SIGF    ???
      COMMON/A26/PIM,PIKPI
      ! A26
      ! ---   PIM     ???
      ! ---   PIKPI   ???
      COMMON/A27/SIGMA,DSIGMA,DELTQ
      ! A27
      ! ---   SIGMA   ???
      ! ---   DSIGMA  ???
      ! ---   DELTQ   ???
      COMMON/A28/BET2,PPI
      ! A28
      ! ---   BET2    ???
      ! ---   PPI     总压
      COMMON/A29/AT2
      ! A29
      ! ---   AT2     ???
      COMMON/A30/C1A,AML1
      ! A30
      ! ---   C1A     ???
      ! ---   AML1    ???
      COMMON/A31/F1,XFK,XFA/A32/ANT,AKG
      ! A31
      ! ---   F1      转子进口流通面积
      ! ---   XFK     ???
      ! ---   XFA     ???
      ! A32
      ! ---   ANT     转速百分比
      ! ---   AKG     空气流量储备系数
      COMMON/A33/QLA1C,HQCK,HQC
      ! A33
      ! ---   QLA1C   ???
      ! ---   HQCK    ???
      ! ---   HQC     ???
      COMMON/A34/PM,P1INIT
      ! A34
      ! ---   PM      第一级入口总压?
      ! ---   P1INIT  压气机进口总压
      COMMON/A35/IJ
      ! A35
      ! ---   IJ      ???
      COMMON/A36/IPEP,KKM,SGM(20),DHQMIN(20),INDEXN
      ! A36
      ! ---   IPEP    ???
      ! ---   KKM     ???
      ! ---   SGM     ???
      ! ---   DHQMIN  ???
      ! ---   INDEXN  ???
      COMMON/A37/WAG1,WAG2,IND8,IND8A
      ! A37
      ! ---   WAG1    ???
      ! ---   WAG2    ???
      ! ---   IND8    ???
      ! ---   IND8A   ???
      COMMON/A38/IND8R
      ! A38
      ! ---   IND8R   ???
      COMMON/A40/C4A
      ! A40
      ! ---   C4A     ???
      COMMON/A41/TTI
      ! A41
      ! ---   TTI     ???
      COMMON/A42/ALF4C
      ! A42
      ! ---   ALF4C   ???
      COMMON/BIS/BIS(30)
      ! BIS
      ! ---   BIS     级间抽气流量或百分比
      COMMON/C1/III,INDEXC,HQKP,HQKPA
      ! C1
      ! ---   III     ???
      ! ---   INDEXC  ???
      ! ---   HQKP    ???
      ! ---   HQKPA   ???
      COMMON/C2/KSI,SL21(6)
      ! C1
      ! ---   KSI     ???
      ! ---   SL21    ???
      COMMON/BLOK/AP1,A1,R,UKR
      ! BLOK
      ! ---   AP1     ???
      ! ---   A1      ???
      ! ---   R       ???
      ! ---   UKR     ???
      COMMON/N/STUP(160),IZY
      ! N
      ! ---   STUP    ???
      ! ---   IZY     ???
      DIMENSION ALF(30),SIG(30),FX(24),W(10)
      DIMENSION STORE(489),GAN(68),BASE(4),PLANE(10,81)
      DIMENSION WCP(10,4),GKI(10),ZC(4)
        INTEGER M
        DOUBLE PRECISION STUP,STORE,GAN,PLANE,PLN1,PLN2,COD
   71 FORMAT(6F13.6)
  102 FORMAT(20I5)
  113 FORMAT(9F12.4)
  116 FORMAT(5X,9HQ(LQMBDA),
     *6X,4HPI K,9X,3HKPD,9X,1HG,10X,4HG PP,6X,6HGPRBIX,6X,3HGPR)
  117 FORMAT(5X,'characteristic axi. stage',5X,'N =',F5.3,
     *5X,'UK =',F5.1,5X,'QL =',F5.3)
  118 FORMAT(5X,5HN OTH,F10.5/2X,4HGKI=,10F10.5)
  119 FORMAT(' 压力不再提高点 ')
      OPEN(25,FILE='PICTURE.s1unknown',ACCESS='DIRECT',RECL=12)
  195 FORMAT(5X,'轮缘结构参数')
  196 FORMAT(5X,10H**********,' 特性线',i2,' 计算',10H**********)
  194 FORMAT(5X,'级号',I2)
  193 FORMAT(8X,'工作叶轮')
  192 FORMAT(8X,'导向器')
  191 FORMAT(8X,'抽气百分比',8X,F13.6)
  190 FORMAT(8X,'计算范围')
  189 FORMAT(8X,'进口参数')
      STEPO=0                 ! 输出step过程
      ALF4C=IWRITE            !记录数组符号
      LSRFIN=0                ! 是否（1/0）寻找与工作状态线的交点
      IF(LSR.EQ.2) THEN
        LSRFIN=1
        LSR=1
      END IF
      IQPZX=IQP
      IF(IQP.EQ.7) IQP=0
      IQP1=IQP
      IQPN=IQP
      IF(IQP.EQ.5) THEN
        IQP=1
      ELSE
        IF(IQP.EQ.6) THEN
          IQPN=5
          IQP=2
        END IF
      END IF
      IZ=IIF      !级数
      N=NNF       ! 前面定义是NNF=1，未知作用（按符号判断可能是特性线数量）
      NTB=NP
      IKBD=IBD
      IF(IKBD.EQ.0) IKBD=1000
      ! 输出轮缘结构参数
      IF(IPE4.EQ.1) THEN
          WRITE(7,196)NOWNV
          WRITE(7,195)
          DO L=IIS,IIF
              WRITE(7,194) L
              WRITE(7,193)
              WRITE(7,71)(A(J,L),J=1,20)      ! 输出各级工作轮中径结构参数
              WRITE(7,192)
              WRITE(7,71)(B(J,L),J=1,20)      ! 输出各级导向器环圈结构参数
              WRITE(7,71)(GKA(J,L),J=1,4)     ! 输出校正系数
              WRITE(7,191)BIS(L)
          END DO
          WRITE(7,190)    ! 打印：计算范围
          DO L=NNS,NNF
              WRITE(7,71)(AN(L,J),J=1,3)      ! 输出各级特性线参数
          END DO
          WRITE(7,71) GB                  ! 输出初始条件？
      END IF
      IPE4=(IPE4+1)/2
      ITT=GB(2)/1000.     !入口总温
      TTI=ITT/10
      GB(2)=GB(2)-ITT*1000
      ! 为了提高精度，转换几何结构单位为mm
      DO L=IIS,IZ
          A(1,L)=A(1,L)*1000.     ! 工作轮入口外缘直径
          A(2,L)=A(2,L)*1000.     ! 工作轮出口外缘直径
          B(13,L)=B(13,L)*1000.   ! 径向间隙
          A(8,L)=A(8,L)*1000.     ! 弦长
          A(13,L)=A(13,L)*1000.   ! 前缘厚度
          A(14,L)=A(14,L)*1000.   ! 前缘到叶栅喉道距离
          B(1,L)=B(1,L)*1000.     ! 导向器入口外缘直径
          B(2,L)=B(2,L)*1000.     ! 导向器出口外缘直径
      END DO
      P1INIT=GB(1)            ! 第一级进口总压
      GB(1)=GB(1)/9.81        ! 第一级入口总压/9.81
      GB(6)=GB(6)/9.81        ! 第一级气体常数
      IF(PB.EQ.0.) PB=1.5
      BASE(1)=IIF-IIS+1       ! 计算总级数
      BASE(3)=A(1,IIS)/1000.  ! 起始级的轮缘直径
      BASE(4)=BASE(3)**2*(1.-A(3,IIS)**2)*0.7854  ! 起始级入口动叶流道面积
      DO L=1,10
          OPTIMS(L)=0.
          BRANCH(L)=0.
          W(L)=GB(L)                  ! 第一行十个数
      END DO
      IF(A(19,IIS).LE.0) THEN         ! 压头特性矫正系数，起始级自动赋值1,后面都为0.9
          A(19,IIS)=1                 ! 压头特性线的旋转系数
          IIS1=IIS+1
          DO L=IIS1,30                ! 后面级旋转系数设为0.9
              A(19,L)=0.9
          END DO
      END IF
      ! 检查修正系数，不允许为0
      DO L=IIS,IIF
          IF(GKA(2,L).EQ.0.) GKA(2,L)=1.
          IF(GKA(3,L).EQ.0.) GKA(3,L)=0.98
          IF(GKA(4,L).EQ.0.) GKA(4,L)=1.
      END DO
      IF(IQP.EQ.2.AND.PR.EQ.0.) IQP=1
      IF(IZAP.NE.0) THEN
          IF(LSR.EQ.0.AND.CU.EQ.0.) CU=20.
          IF(LSR.NE.0.AND.ABS(WCP(1,1)).LE.0.) THEN
              LSR=0
              IF(CU.EQ.0.) CU=20.
          END IF
      END IF
      IF(IQP+IZAP.GT.0.AND.EP.EQ.0.) EP=0.05    ! ？
      IF(ITT.EQ.0) TTI=GB(2)
      PA=57.296                                 !180/3.14
      MSC=0                                     !?
      IQPS=IQP                                  
      IZ1=0.25*IIF
      IZ2=0.5*IIF
      IZ3=0.75*IIF
      IXIG=0                                   !与流量极值有关
      QQS1=QQS                                 !最左侧流量系数
      DO 6 M=NNS,NNF      ! 每个特性线的计算,索引为M
      IQPZZ=IQP
      IF(IXIG.EQ.1) THEN
          QQS=0
          QQF=0
      ENDIF
      IXIG=0
 6475 CONTINUE
      STEQ1=0.
      STEQ2=100.
      ISTEQ=0
      AKPDMA=0.
      PMIN=100.
      FAN=AN(M,2)
      PP=0.
      AKPX=0.
      SDELTQ=0.
      QH=0.
      QGP=0.
      PGP=0.
      NG=0
      IG=1
      IGG=0
      LAR=6+(IIS-1)*2   !控制级循环
      INDEX1=0
      IQP=IQPS
      IQ=0
      DO L=1,10
          GKI(L)=GKA(1,L)   ! 落后角修正系数，GKI不同于其他3个修正系数
          GB(L)=W(L)
      END DO
      INDEXN=0
      INDEXC=0
      INDEXQ=0
      IND8A=0
      IND8R=0
      IPEP=0
      IGG=0
      ILR=0
      MW=0
      MV=0
      PW=0
      IDQC=0
      POINTS=0.
      INDEXB=0                !
      INDEXD=0                !
      INDEXZ=0                !
      DSIGMA=DSIG             !垂直段总压系数变化步长
      IF(DSIG.LT.1E-6) DSIGMA=0.015
      UKP=GB(3)*AN(M,1)       ! 第一个工作轮前缘圆周速度 * 特性线转速百分比
      UKM=UKP
      KG=0
      IPES=0
      IF(IQP.EQ.3) INDEXB=4
      DO L=1,40
          KRAN(L)=KPAH(L)     ! 控制是否输出数组数据
      END DO
      DO L=1,81
          PLANE(M,L)=0.       ! 填充0到这个数组
      END DO
      INDEX2=0
      ILMN=0
      INDEX4=0
      INDEX6=0
      INDEX7=0
      INDEX8=0
      INDEX9=0
      GUIDE=0
      IL=IIS
      IM=0
      IN=0        ! 啥意思啊
      KM=0
      KMI=0
      K=0
      PIKPI=0
      PIM=0
      III=0
      JJJ=0
      KSI=0
      IF(QQS.GT.1E-6) AN(M,2)=QQS         ! 流量计算初始点
      IF(AN(M,1).LE.0.78) GB(8)=2*GB(8)   ! 转速百分比小于0.78时流量计算步长为2倍
      DELTQ=GB(8)       !特性计算步长，沿q
   57 CONTINUE
      J=10
      IF(K.EQ.1) J=8
      DO L=1,J
          GB(L)=W(L)
      END DO
      FX(22)=1.
      STAGE(22)=1.
      TT1=GB(2)       ! 第一级入口处总温
      PP1=GB(1)       ! 第一级入口处总压
      ARFA1=GB(4)     ! 第一个工作轮入口绝对气流角
      TTT=TT1
      PPP=PP1
      GO TO 10
   56 CONTINUE
      IF(INDEX6.EQ.0) GOTO 57
      IF(MV.EQ.1) GOTO 2101
      IF(INDEX9.NE.1) THEN
          IN=GUIDE/10
          DSIGMA=0.01
          IF(IN.GT.IZ1) DSIGMA=0.015
          IF(IN.GT.IZ2) DSIGMA=0.02
          IF(IN.GT.IZ3) DSIGMA=0.025
          IF(DSIG.GT.0.) DSIGMA=DSIG
      END IF
      CALL WAY2(STAGE,FX,GB)          ! 函数功能未知
      IF(KMI.EQ.2) KSI=0
      STAGE(22)=PPI
      FX(22)=PPI
      !                                   /************ 级循环，行号 1 ************/
  10  DO 1 I=IIS,IIF
      INDEXC=0
      MG=IIS/I
      IDQC=(1-MG)*IDQC
      CALL PE4ATI(IPE4,KPAH,FKPAH,KRAN,I,AN,M,PP,AKPX,IG)     ! 检查参数范围是否符合要求，并生成输出控制数组
      UKZ=UKP
      IF(AN(M,3).NE.0.) THEN
          TKT=1
          UKZ=UKM*MG+UKP*(1-MG)       ! ?????
          IF(I.EQ.IKBD) THEN
              IF(NTB.NE.0) TKT=SQRT(GB(2)/TT1)    ! 此处NTB=0
              UKZ=UKP*AN(M,3)*TKT
          END IF
      END IF
      AP1=GB(1)           ! 入口总压
      A1=GB(4)            !第一个工作轮入口绝对气流角
      II=I
      IJ=0
      IZW=IIS*8-7             ! 
      IZY=IZW*MG+IZY*(1-MG)   ! 
      ! 获取等熵过程指数
      IF(GB(7).GT.0.) THEN
          AK=GB(7)            
      ELSE
          AT1=GB(2)           ! 工作轮进口总温
          AK=TEMP(AT1)        !求等熵过程指数
      END IF
      ILMN=IN-II
      IF(ILMN.GT.0) GO TO 816     ! 跳过第IN级
      GAU=GB(1)*9.81
      GAR=GB(4)           ! 在第一级工作轮进口处的气流角
      GAN(LAR)=COD(GAU,0,GB(2),2,5)
      UK1=UKZ*A(1,II)/A(1,IIS)    ! 工作轮前缘圆周相对速度
      FORM1=ALF(I)        ! 出口气流角
      FORM2=SIG(I)        ! 级后总压恢复系数
      DO K=1,20
          RK(K)=A(K,I)
          HA(K)=B(K,I)
      END DO
      HA(17)=B(17,1)      ! 是否计算雷诺数
      RK(5)=RK(5)/PA      ! 角度值换算为弧度值
      RK(6)=RK(6)/PA
      HA(5)=HA(5)/PA
      HA(6)=HA(6)/PA
      RK(16)=RK(16)/PA
      HA(16)=HA(16)/PA
      ! 速度足够小的时候忽略叶型尖角的影响
      IF(UKM.LE.150.AND.RK(16).NE.0.) THEN
          RK(10)=0
          RK(16)=0.
      END IF
      IF((INDEX4.LE.0).OR.(I.GT.IIS)) THEN
          CALL OPTIMA(II,RK,HA,RKI,HAI,
     *GB,GKA,X,Y,GDK,GDA,STAGE,INDEX1,INDEX4,OPTIMS,GKI,M)      ! 计算能头、效率等
          IF(INDEX4.LE.0) THEN
              DO L=1,10
                  ARKI(II,L)=RKI(L)
                  AHAI(II,L)=HAI(L)
              END DO
              AHAI(II,11)=STAGE(19)
          END IF
          GB(4)=GB(4)*PA      ! 弧度值换算为角度值
      ELSE
          DO L=1,24
              STAGE(L)=FX(L)
              C(1+L/7)=ZC(1+L/7)
          END DO
          DSITAK=W1
          DSITAA=W2
          DO L=1,10
              RKI(L)=ARKI(II,L)
              HAI(L)=AHAI(II,L)
          END DO
      END IF
      CALL NEOPT (M,II,INDEX2,GB,AN,
     *RK,HA,RKI,HAI,STAGE,GKA,GDK,GDA,FFX)
      IF(INDEXC.GT.0.AND.A(3,II).GT.0.73) IDQC=II
      GAN(LAR+1)=COD(GAR,3,C1A,4,6)       ! GAN
      LAR=LAR+2
      XQLA1C=QLA1C
      XHQCK=HQCK
      XHQC=HQC
      QLA1C=0
      HQCK=0
      HQC=0
      IF(I.EQ.IIS) THEN
          ! 压气机第一级的参数处理
          AAAA=ALFA1
          QL1=AN(M,2)
          D=SQRT(AK*(2./(AK+1.))**((AK+1.)/(AK-1.)))
          D=D*SQRT(9.81/RRR)
          F1=3.1415*DK1**2*(1.-DDK1**2)*0.000001/4.
          G=D*F1*PPP*QL1*SIN(AAAA)/SQRT(TTT)          ! 流量计算公式
          IF(AKGK.GT.1E-2) G=G*AKGK
          GPP=G*10328.746/PPP*SQRT(TTT/288.)
      END IF
      IF(FORM1.GT.1E-6) GB(4)=FORM1
      IF(ALF(1).GT.1E5) GB(4)=HAI(2)
      ! 级后总压恢复系数对气体总压和压比的影响
      IF(FORM2.GT.0.000001) THEN
          GB(1)=GB(1)*FORM2
          STAGE(22)=STAGE(22)*FORM2
      ENDIF
      IF(INDEX4.LT.1.AND.I.EQ.IIS) THEN
          DO L=1,24
              FX(L)=STAGE(L)
              ZC(1+L/7)=C(1+L/7)
          END DO
          W1=DSITAK
          W2=DSITAA
      ENDIF
      GB(4)=GB(4)*PA
      IF(IPES.EQ.II) THEN
          INDEX2=1
          INDEX8=5
          IND8=0
          IND8A=1
          IPES=0
          AN(M,2)=AN(M,2)+ABS(DELTQ)                  ! 流量循环？
      ENDIF
      IF(INDEX2.GT.0) GOTO 805
      STUP(IZY+5)=COD(STAGE(15)*PA,2,STAGE(2)*PA,2,6)
      GOTO 1
  816 LAR=LAR+2                                       ! LAR参数：
    1 IZY=IZY+8
      !                                       /************ 级循环结束，行号 1 ************/
      IF(IN.LT.99) GOTO 2101
      MV=1
      PPI=PMIN
 2101 IZY=8*IIS-7
      NG=NG+1
      WAG2=WAG1
      PPI=STAGE(22)*(1-MV)+(PPI-0.2)*MV
      GB(2)=GB(2)*(1-MV)+TW*MV
      IF(PMIN.GT.PPI) PMIN=PPI
      BN=AN(M,1)
      QL1=AN(M,2)
      PP=PPI
      CN=AN(M,1)*GB(3)*SQRT(288/TT1)                  ! 换算转速？
      STEQ1=PPI
      IF(ABS(STEQ1-STEQ2).LT.1E-5) THEN
          ISTEQ=ISTEQ+1
          IF(IPE4.NE.0) THEN
              ! 可能有一些输出参数
          END IF
      END IF
      STEQ2=PPI
      IF(ISTEQ.EQ.5) THEN
          IF(INDEX6.GT.0) GOTO 66
          ISTEQ=1
          DELTQ=0.2
      END IF
      IF(GB(7).EQ.0) THEN
          AK=TEMP((GB(2)+TT1)/2.)
          CALL AKPDI(TT1,GB(2),PP,GB(6)*9.81,AKPX)       !!!!!!!!!!!! 计算等熵效率
      ELSE
          AKPX=(TT1*(PP**((AK-1.)/AK)-1.))/(GB(2)-TT1)   !!!!!!!!!!!! 等熵效率
      END IF
      ! 计算出口物理流量
      QQQ=QLA1*F1*1E-6*D/SQRT(9.81/RRR)/SQRT(RRR*9.81)
      QQQ=QQQ*GB(1)*9.81/SQRT(GB(2))*SIN(ALFA1)
      !QQQ=QLA1*F1/1000**2*AM/SQRT(RRR*9.81)*GB(1)*9.81/
      !*PI/SQRT(AT1)/AKG*SIN(ALFA1)
      QK=0
      DO LI=1,IIF
          QK=QK+BIS(LI)
      END DO
      !QQQ=G*(1.-QK)
      GPRBIX=G*(1.-QK)*10328.746/GB(1)*SQRT(GB(2)/288.)           ! 轴流压气机出口折合流量
      IF(AKPDMA.LE.AKPX) THEN
          AKPDMA=AKPX
          STORE(327)=GPP
          STORE(328)=PP
          STORE(329)=AKPX
          LJ=0
          DO LI=330,489        ! STORE的第330至489号用于存放STUP数组
              LJ=LJ+1
              STORE(LI)=STUP(LJ)
          END DO
      ENDIF
      ! 是否打印一些输出信息
      IF(IPE4.NE.0) THEN
          IF(NG.EQ.1) then
              WRITE(7,117)BN,CN,QL1
              write(numberchar,'(f5.3)') bn
              if(outflag.eq.1) then
              write(7002,*)'zone t="character'//trim(numberchar)//'"'
              endif
          endif
          IF(KRAN(39).NE.0.OR.KRAN(40).NE.0) THEN
              WRITE(7,116)
              WRITE(7,113)QL1,PP,AKPX,G,GPP,GPRBIX,QQQ
          ELSE
              IF(NG-1.EQ.((NG-1)/10)*10) WRITE(7,116)
              WRITE(7,113)QL1,PP,AKPX,G,GPP,GPRBIX,QQQ
              if(outflag.eq.1) then
                  WRITE(7002,*) G,PP,AKPX
              endif
              NPOINT=IG
              AGP(1)=GPP
              AGP(2)=PP
              AGP(3)=AKPX
              WRITE(25,REC=NPOINT) AGP
          ENDIF
      ENDIF
      IF(IIF.EQ.IIS.AND.IXIG.EQ.0.AND.IQP.EQ.0) THEN  ! 单级压气机
          IF(XHQC.GT.50.) THEN
              AN(M,2)=AN(M,2)-0.05
          ELSE
              QQS=XQLA1C
              QQF=1.
              AN(M,2)=XQLA1C
              IQPZZ=7
              IXIG=1
              WRITE(7,*) ' QLA1C= ',XQLA1C,' HQC= ',XHQC,' HQCK= ',XHQCK
              outflag=1
          ENDIF
          GOTO 6475
      ENDIF
      IF(PP.LE.0.5) GO TO 66
      IF(MV*IQP.GT.0) GOTO 66
      IF(INDEX6.EQ.0) IGG=IG
      IF(IG.NE.IGG.AND.IGG.NE.1) THEN
          NLK=IG-2
  845     DO KU=1,NLK
              NKU=2*KU
              NKU2=NKU+2
              CALL DECOD(PLANE(M,NKU),5,RAS,4,QAS,5)
              CALL DECOD(PLANE(M,NKU2),5,RAF,4,QAF,5)
              IF(RAS.GT.RAF) THEN
                  PLN1=PLANE(M,NKU)
                  PLN2=PLANE(M,NKU+1)
                  PLANE(M,NKU)=PLANE(M,NKU2)
                  PLANE(M,NKU+1)=PLANE(M,NKU2+1)
                  PLANE(M,NKU2)=PLN1
                  PLANE(M,NKU2+1)=PLN2
              ENDIF
          ENDDO
          NLK=NLK-1
          IF(NLK.GT.1) GO TO 845
          IF(IQPZZ.NE.7) THEN
              PGR=0
              DO KU=1,IG
                  NKU=2*KU+1
                  CALL DECOD(PLANE(M,NKU),5,PAS,5,CAS,5)
                  IF(PAS.GT.PGR) THEN
                      PGR=PAS
                      JLK=KU
                  ENDIF
              ENDDO
              IGJ=IG-JLK+1
              DO KU=1,IGJ
                  NKU=2*KU
                  NKUJ=(KU+(JLK-1))*2
                  PLANE(M,NKU)=PLANE(M,NKUJ)
                  PLANE(M,NKU+1)=PLANE(M,NKUJ+1)
              ENDDO
              IG=IGJ
          END IF
      ENDIF
      GAN(1)=AN(M,1)  ! 转速百分比
      GAN(2)=GPP      ! 空气折合流量
      GAN(3)=QL1      ! 流量系数
      GAN(4)=PP       ! 压气机增压比
      GAN(5)=AKPX     ! 等熵效率
      GAU=GB(1)*9.81  ! 总压
      KG=KG+1
      GAN(LAR)=COD(GAU,0,GB(2),2,5)
      GAN(LAR+1)=COD(GB(4),3,C4A,4,6)
      LAR=6
      IF(AKPX.LT.0.) AKPX=0.001
      PLANE(M,2*IG)=COD(GPP,4,QL1,5,5)
      PLANE(M,2*IG+1)=COD(PP,5,AKPX,5,5)
      BIG=IG
      PLANE(M,1)=COD(AN(M,1),5,BIG,0,2)
      IG=IG+1
      IF(INDEXZ.EQ.1) IG=IG-1
      IF(IG.GT.40) IG=40
      IF(QQF.GT.1E-6.AND.AN(M,2).GE.QQF) GO TO 66
      IF(INDEXN.NE.0.OR.MV.NE.0) THEN
          IF(MV.NE.1) THEN
              L=GUIDE/10+IIS
              IF(INDEX6.EQ.1) L=GUIDE/10+1
              IF(JJJ.GT.0) L=JJJ+1
              IF(L.GT.IZ) L=IZ
              DHX=DHQMIN(L)
              IPEP=L
              JJJ=0
              DO LP=L,IZ
                  IF(DHX.GT.DHQMIN(LP)) THEN
                      DHX=DHQMIN(LP)
                      IPEP=LP
                  END IF
              END DO
              IF(DHX.GT.0.) MW=1
          END IF
          PW=PP*MW+PW*(1-MW)
          TW=GB(2)
          IF(PW/PP.GE.PB) GOTO 66
          IM=IPEP*(1-MW)
          KM=SGM(IPEP)*(1-MW)
          MV=MW+MV
          MW=0
          INDEXN=0
          IF(INDEX6.EQ.0) IPES=IPEP
          IF(IPES.EQ.1) IPEP=0
          K=1
          GOTO 56
      END IF
      GB(1)=PP1
      GB(2)=TT1
      GB(4)=ARFA1
      UKP=UKM
      ! STEP 过程
      IF(IQP.LE.0.OR.IQP.GT.3) THEN
          CALL STEPIN(INDEXZ,INDEX2,INDEX6,INDEX9,INDEXB,LSR,IM,KM,IQP,
     *PR,IZAP,CU,ZAP,AN(M,2),GPP,PP,DELTQ,SIGMA,DSIGMA,QSET,DQSET,SISET,
     *DSISET,SDELTQ,WCP,EP,GUIDE,SGUIDE,KSTEP)
          IF(STEPO.EQ.1)WRITE(7,*)'STEPIN'
          IF(IQP.EQ.3.OR.INDEXZ.EQ.2) THEN
              ! INDEXZ=1
              IF(INDEXZ.EQ.2) THEN
                  MG=INDEXZ/2
                  INDEXZ=1+2*MG
              END IF
              STORE(1)=GPP
              STORE(2)=PP
              STORE(3)=AKPX
              DO LK=4,163
                  STORE(LK)=STUP(LK-3)
              END DO
          END IF
          MG=(QGP+4)/(AN(M,2)+4)*INDEXB/4
          AN(M,2)=AN(M,2)+(QGP-AN(M,2))*MG
          QL1=AN(M,2)
          IF(KSTEP.EQ.3) STAGE(22)=1.
          IF(LSRFIN.EQ.1.AND.INDEXZ.EQ.3) KSTEP=4
          ! KSTEP 步骤说明
          ! ---   1:  820 STEPAT
          ! ---   2:  56  重新开始级循环，包括WAY2
          ! ---   3:  10  重新开始级循环
          ! ---   4:  66 特性线计算循环结束
          GO TO (820,56,10,66),KSTEP
      END IF
 820   CALL STEPAT(INDEX2,INDEXB,INDEXD,INDEXQ,INDEX6,INDEX8,
     *IND8R,QQS,QQF,GB(8),ES,AN(M,2),PP,GPP,Q1,P1,QGP,PGP,GGP,Q2,P2,
     *DELTQ,QMEM,DQMEM,ZAP,IQP,KSTEP)
      IF(STEPO.EQ.1)WRITE(7,*)'STEPAT'
      IF(INDEXB*IZAP.EQ.3) QH=QGP
      IF(IQP.EQ.1.OR.PP.GE.PGP) THEN
          STORE(164)=GPP
          STORE(165)=PP
          STORE(166)=AKPX
          DO KL=167,326
              LK=KL-166
              STORE(KL)=STUP(LK)
          END DO
      END IF
      IF(INDEXB.EQ.3) THEN
          IF(LSR.NE.0) THEN
              IF(PGP.LE.GWCP(GPP,WCP)) ILR=1
          END IF
          QL1=AN(M,2)
      END IF
      INDEXZ=INDEXZ*(1-ILR)+2*ILR
      ! KSTEP 步骤说明
      ! ---   1:  805 STEPON
      ! ---   2:  10  重新开始级循环
      GO TO (805,10),KSTEP
  805 CALL STEPON(M,IQP,IQ,IDQC,INDEXQ,INDEX2,INDEX4,INDEX6,INDEX7,
     *IPE4,GB,AN,DELTQ,PP,QL1,AKPX,QH,QGP,PGP,QQS,QQF,IM,
     *KM,INDEX8,IND8R,INDEXB,KSTEP,BRANCH,POINTS)
      IF(STEPO.EQ.1)WRITE(7,*)'STEPON'
      MG=1-KSTEP/5
      IG=MG*(IG-1)+1
      IGG=IGG*MG
      LAR=6
      ! KSTEP 步骤说明
      ! ---   1:  9   下一步并重新开始级循环
      ! ---   2:  232 跳过IQP判断
      ! ---   3:  66 特性线计算循环结束
      ! ---   4:  66 特性线计算循环结束
      ! ---   5:  57  重新开始级循环，包括总温总压赋值及WAY2部分
      GO TO (9,232,66,66,57),KSTEP
    9 IF(IQP.GE.2) INDEX1=1
          FX(22)=1
      GO TO 10
  232 CONTINUE
      IF(IJ.GT.0.AND.JJJ.EQ.0) JJJ=III
      CALL WAY1(M,INDEX2,GB,AN,PPI,SIGMA,DELTQ,DSIGMA,K,FX,FFX,PB)
      LAR=6
      MG=2/(IM+1)
      MG=ABS(2-MG)/2
      JJJ=JJJ*MG
      IF(IQP.LT.2) GOTO 861
      IF(ABS(DSIGMA).LE.1.E-6) GOTO 888
      IF(INDEX8.NE.2) GO TO 861
      INDEX8=4
      AN(M,2)=AN(M,2)-2.*DELTQ
  888 IF(PP+EP.GE.PR) GOTO 232
      IF(IPE4.GT.0) WRITE(7,119)
      POINTS=2
      GO TO 66
  861 CONTINUE
      INDEX2=0
      IF(K.EQ.1) GO TO 56
      MSC=MSC+1
   66 CONTINUE
      IF(OPTIMS(M).GT.0.AND.OPTIMS(M).LT.100.) WRITE(7,118) AN(M,1),GKI
      IF(OPTIMS(M).GT.0.AND.OPTIMS(M).LT.100.) WRITE(0,118)AN(M,1),GKI
      IF(ILR.EQ.1) STORE(1)=999999
      IF(ILR.EQ.1) WRITE(7,198)
      IF(ILR.EQ.1) WRITE(0,198)
      AN(M,2)=FAN
      DO L=1,10
          GB(L)=W(L)
      END DO
      BASE(2)=MSC
    6 CONTINUE                ! 特性线计算结束
      QQS=QQS1
      DO L=IIS,IIF
          A(1,L)=A(1,L)/1000.
          A(2,L)=A(2,L)/1000.
          A(8,L)=A(8,L)/1000.
          A(13,L)=A(13,L)/1000.
          A(14,L)=A(14,L)/1000.
          B(1,L)=B(1,L)/1000.
          B(13,L)=B(13,L)/1000.
          B(2,L)=B(2,L)/1000.
      END DO
      DO L=1,10
          GB(L)=W(L)
      END DO
      GB(1)=GB(1)*9.81
      GB(6)=GB(6)*9.81
      IQP=IQPZX
      GB(2)=GB(2)+ITT*1000
  198 FORMAT('line 1559: unknown label text 198')
      RETURN
      END