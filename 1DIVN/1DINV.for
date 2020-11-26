
      PROGRAM OCK

      !数组定义，详见 notion 说明
      DIMENSION KPAH(40),FKPAH(10),LCP(10,4)
     *,GAN(68),STORE(489),BASE(4),GKI(10),PLANE(10,81),AA(30,65)
      DIMENSION X(16),Y(16),D3(9),D31(30),D32(30),ID1(33),X2(16),Y2(16),
     *ID2(24),GKA(4,30),AN(10,3),GB(10),A(20,30),B(20,30),X3(16),Y3(16)
     *,AN2(10,3)
      DIMENSION SIG(30),ALF(30),BRANCH(10),OPTIMS(10)
      DIMENSION DI(30,2),SI(30),CYMP(20),KP(30,3)
      DIMENSION TETA(10,30)
      DIMENSION BXD1(10),BXD2(50)
      REAL ID1,LCP,ID2
      DOUBLE PRECISION  STORE,GAN,PLANE
      DIMENSION U(10),GI(20),PII(20),ETI(20),PLA(4,40),PIOC(10,40)
     *,GPROC(10,40),ETOC(10,40),SRN(3,30),PIC(40),ETC(40)
      DATA X/.0,.02,.05,.1,.15,.2,.25,.3,.4,.5,.6,.7,.8,.9,.95,1./
     *,Y/.0,.0175,.02625,.03475,.03975,.04425,.04625,.04875,.05,.04875,
     *.04525,.038,.028,.01675,.00925,.0/                                ! 叶型厚度分布
      DATA X2/.0,.0125,.025,.05,.075,.1,.2,.3,.4,.5,.6,.7,.8,.9,.95,1./
     *,Y2/.0,.0165,.0227,.0308,.0362,.0402,.0483,.05,.0489,.0457,.0405, 
     *.0337,.0254,.016,.0106,.0/                                        ! 叶型厚度分布                    
      DATA X3/.0,.0125,.025,.05,.075,.1,.2,.3,.4,.5,.6,.7,.8,.9,.95,1./
     *,Y3/.0,.01169,.01574,.02177,.02647,.0304,.04143,.0476,.04996,
     *.04812,.04146,.03156,.01987,.0081,.00306,.0/                      ! 叶型厚度分布


      OPEN(7,FILE='Result.1D')        ! 输出文件，文件号为 7
C     open(7006,FIlE='parastos2.s12s2')
	OPEN(7006,FILE='Par_AVR.1DtoS2')
      OPEN(4,FILE='ax18.dat')     ! 输入原始数据文件

      IR=2                            ! 设计计算控制参数（1计算；2设计）
      IF(IR.EQ.1) GOTO 56             ! 进行计算
      ! 设计方案
      CALL BBOD2(IOCK,PIK,GPRR,GB,GKA,AN,A,B,IZ,KF,NV,TETA,
     *OBC,ETCBK,D2C,BET2L,ALA3)
      SVCBK=1.
      BXD1(4)=BET2L/57.296
      BXD1(1)=PIK
      G3=GPRR
      IF(IR.EQ.2) GOTO 55             ! 计算方案跳转至55
      IF(IR.NE.1.AND.IR.NE.2) GOTO 33
   56 CONTINUE                        ! 开始计算结构参数
      CALL BBOD1(OB,PK,G3,PBX,T1,XTA3,
     *AKP,IZ,K,KF,KDB,UBX,PBHA,
     *D3,D31,ID1,IZ2,BXD1,ALF1,IOCK,NV,AN,TETA,GB,PSI)    ! 读入原始数据，并打印计算结果输出列表的第一大块
      IF(OB.LT.900.) OB=60.*OB/3.14159/D31(1)             ! 处理转速数据
      OPR=OB*SQRT(288./T1)                                ! 折合转速
      ! TODO: 可以优化此处的GOTO结构（合并101和102）
      IF(IOCK.EQ.2) GOTO 102                              ! 运算方案符号2--离心级计算，跳转至102处
      CALL PCPM(OB,PK,G3,PBX,T1,XTA3,AKP,IZ,K,KF,KDB,
     *UBX,PBHA,D3,D31,D32,ID1,ID2,A,B,GKA,GB,KP,SI,DI,
     *AA,CYMP,IZ2)                                        ! 
      GPRR=G3/PBHA/UBX*SQRT(T1/288.)
  102 CONTINUE
      IF(IOCK.EQ.2) GOTO 101
      CALL ARCH1(SI,DI,AA,IZ,CYMP,OB)
  101 CONTINUE
      ALF1=ALF1/57.296
      IF(IOCK.EQ.2.OR.IOCK.EQ.1) GOTO 100
      CALL PERKSI(AA,BXD1,IZ,PKSI,DLP,RV,DLC,GAM,ALP,ALF1)
  100 CONTINUE
      IF(IOCK.EQ.1) BXD1(1)=1.
      IF(IOCK.EQ.1) GOTO 22
      ! 开始计算离心压气机部分
      PIK=BXD1(1)
      SVCBK=1.
      A1C=BXD1(8)/57.296
      IF(ABS(A1C-ALF1).GT..0001) SVCBK=.995
      IF(IOCK.EQ.2) PK=1.
      IF(IOCK.EQ.2) PKSI=1.
      IF(IOCK.EQ.2) AA(IZ+1,26)=T1
      BXD1(1)=BXD1(1)/PK/PKSI/SVCBK
      POC=PBX*PK*PKSI*SVCBK
      BXD1(4)=BXD1(4)/57.296
      OBC=OB*BXD1(9)
      G3C=G3*BXD1(10)
      CALL CENTRO(BXD1(4),A1C,BXD1(7),BXD1(6),OBC,G3C,AA(IZ+1,26),POC,
     *BXD1(1),IERR,D2C,H2C,ALRKC,DGC,PCBK,TCBK,ETCBK,H3C,D3C,H4C,D4C,
     *D5VC,ALA3,PSI,GAMLD)
      CYMP(1)=PK
      IF(IOCK.EQ.12) ETAS=(PIK**.2857-1.)/((1.+(PK**.2857-1.)/
     *CYMP(2))*(1.+(BXD1(1)**.2857-1.)/ETCBK)-1.)
      IF(IOCK.EQ.2) ETAS=0.
      IF(IOCK.EQ.1) ETAS=0.
      BXD2(35)=H4C
      BXD2(26)=D3C
      BXD2(16)=D2C
      BXD2(17)=H2C
      BXD2(30)=H3C
      BXD2(34)=D4C
      BXD2(41)=D5VC
      BXD2(40)=DGC
      BXD2(6)=ALRKC
   22 CONTINUE
      ! 离心压气机部分计算完成
      CYMP(1)=PK
      CALL PER(AA,BXD1,BXD2,IZ,DLP,G3,RV,PIK,GAM,ETAS,SI,DI,CYMP,
     *PKSI,ALP,ETCBK,A1C,ALF1,IOCK,OPR,PSI,GAMLD)
      IPR=2                               ! 进行特性计算
      IF(IPR.EQ.1) GOTO 33                ! 特性计算控制参数（1为不进行计算，结束程序）
      IF(IOCK.EQ.2) GOTO 31               ! 离心压气机不计算特性
      ! 计算结构参数
      DO 3555 I=1,IZ
      A(20,I)=0.                                                  ! 在分离边界上，工作轮空气动力载荷的判据值
      A(13,I)=.16*A(8,I)*A(12,I)*A(10,I)                          ! 型面前缘的厚度
      A(16,I)=ATAN(1.2*A(12,I))*57.3*A(10,I)                      ! 型面尖角的一半
      A(14,I)=A(8,I)/A(7,I)*COS((A(5,I)-A(16,I))/57.3)*A(10,I)    ! 由型面前缘到叶栅喉道的距离
      B(17,I)=0.                                                  ! 是否计算雷诺数（1/0）
      GKA(1,I)=1.                                                 ! 落后角校正系数
      GKA(2,I)=1.                                                 ! 空气流量储存系数
      GKA(3,I)=.98                                                ! 理论压头减少系数
      GKA(4,I)=1.                                                 ! 最大效率减少系数
 3555 CONTINUE
   55 CONTINUE
      ! 开始计算气动参数
      IF(IOCK.EQ.2) GOTO 31
      ! 将亚音速级初始叶片厚度分布传入X,Y两个数组，初始为 BC-10叶形
      IF(KF.EQ.3) THEN
          DO I=1,16                       ! NACA-65叶形
          X(I)=X3(I)
          Y(I)=Y3(I)
          END DO
      ELSEIF (KF.EQ.2) THEN
          DO I=1,16                       ! C-4叶形
          X(I)=X2(I)
          Y(I)=Y2(I)
          END DO
      END IF
      ! 
      GB(7)=1.4                           ! 等熵指数K=1.4
      ALF(1)=1000000.                     ! 第一级级后气流角度
      AKGK=1.                             ! 整个压气机流量系数
      DO 3556 I=1,40                      ! 默认不输出所有级的参数???
      KPAH(I)=0
 3556 CONTINUE
      DO 3560 I=1,30                      ! 级后总压恢复系数设为0
      SIG(I)=0
 3560 CONTINUE
      DO 3561 I=1,10                      ! 默认不输出任何参数????
      FKPAH(I)=0
 3561 CONTINUE
      DO 3562 I=1,4                       ! TODO: 找出LCP数组的作用
      LCP(10,I)=0
 3562 CONTINUE
      ! 应该是一些参数的初始化
      NTB=0       !
      IKBD=0      ! 可能是前几级的某个参数计算阈值
      IQP=0
      IWRITE=0
      NNS=1
      NNF=1
      IIS=1
      IIF=IZ
      QQS=0.
      QQF=0.
      IPE4=1      ! 控制是否输出参数及说明
      DSIG=0.
      IZAP=0
      LSR=0
      CU=0.
      EP=0.
      PR=0.
      PB=1.5
      ! TODO: 搞清楚上面的代码是用来做什么的
   31 CONTINUE
      OPEN(7002,FILE='character_mepic.tec')
      write(7002,*) 'VARIABLES = "Q","p_ratio","efficiency"'
      ! 接下来到3557是一个大的循环，计算每一条特性线，NV为特性线数量
      DO 3557 I=1,NV
c      write(7002,*)'zone t="character'//'"'
      AN2(1,1)=AN(I,1)
      IF(IOCK.EQ.2) GOTO 32       ! 离心压气机不用计算以下部分
      ! 至3558处计算气流角
      GB(4)=GB(4)+TETA(I,1)
      AN2(1,2)=AN(I,2)
      AN2(1,3)=AN(I,3)
      DO 3558 J=1,IZ              ! 计算导向器叶片出入口处的结构角
      B(5,J)=B(5,J)+TETA(I,J+1)
      B(6,J)=B(6,J)+TETA(I,J+1)
 3558 CONTINUE
      ! 调用HARIKA算法，对应说明书中子程序AXIS_COMP_PROPER
       CALL HARIKA(GB,GKA,AN2,A,B,X,Y,SIG,ALF,NTB,IKBD,
     *IQP,IWRITE,AKGK,NNS,NNF,IIS,IIF,QQS,QQF,IPE4,KPAH,FKPAH,
     *DSIG,IZAP,LSR,LCP,CU,EP,PR,PB,BRANCH,OPTIMS,GKI,
     *POINTS,BASE,STORE,GAN,PLANE,I)
      WRITE(0,196)I
  196 FORMAT(2X,'I= ',I2)
      CALL DECOD(PLANE(1,1),2,A3,5,B3,0)
      LB3=B3
      DO 8 L=1,LB3
      CALL DECOD(PLANE(1,2*L),5,A1,4,B1,5)
      CALL DECOD(PLANE(1,2*L+1),5,A2,5,B2,5)
      PLA(1,L)=A1
      PLA(2,L)=A2
      PLA(3,L)=B2
      PLA(4,L)=B1
      IF(A2.LT.1.) PLA(3,L)=B2-1.
    8 CONTINUE
      DO 6 IS=1,IZ
      CALL DECOD(STORE(167+8*(IS-1)),5,A4,4,B4,4)
      CALL DECOD(STORE(168+8*(IS-1)),6,A5,4,B5,4)
      CALL DECOD(STORE(171+8*(IS-1)),6,A6,4,B6,3)
      SRN(1,IS)=A4
      SRN(2,IS)=A6
      SRN(3,IS)=A5
    6 CONTINUE
      IF(IOCK.EQ.1) GOTO 3
      ! 轴流压气机不用计算以下至3部分
      RN=OBC*SQRT(288./AA(IZ+1,26))
      GR=G3C/PK/PBHA/UBX/PKSI*SQRT(AA(IZ+1,26)/288.)
      GOTO 34
      ! 离心压气机不用计算以上部分
   32 CONTINUE
      RN=OBC          ! 离心轮转速与轴流转子转速之比
      GR=G3/SVCBK
      LB3=1
      ! 猜测RN和GR是有关离心组件的参数
   34 CONTINUE
      PIR=BXD1(1)
      ETR=ETCBK
      R=287.
      TB=288.
      AK=1.4
      D2=D2C
      BT2=BXD1(4)*57.296
      K1=1
      K2=1
      K3=1
      KK=1
      NOM=1
      IF(IOCK.EQ.12) WRITE(7,117)A3
      IF(IOCK.EQ.12) WRITE(7,116)
      IIQP=0
      IDOP=0
      ITC=0
      DO 2 M=1,LB3
      PIC(M)=0.
      ETC(M)=0.
      IF(IOCK.EQ.2) GOTO 35
   15 CONTINUE
      TOV=288.*((PLA(2,M)**.2857-1.)/PLA(3,M)+1.)     ! 离心轴流参数
      GPRV=PLA(1,M)/PLA(2,M)*SQRT(TOV/288.)           ! 离心轴流参数
      TNZ=OBC*AN2(1,1)*SQRT(288./TOV)                 ! 离心轴流参数
   35 IF(IOCK.EQ.2) TNZ=OBC*AN2(1,1)
      IF(TNZ.LE.RN) NOM=0
      IF(TNZ.GT.RN) KK=0
      CALL HARCK(TNZ,GR,PIR,ETR,RN,ALA3,TB,AK,D2,BT2,K1,K2,K3,
     *R,U,KK,GI,PII,ETI,NOM)
      IF(IOCK.EQ.2) GOTO 36
      ! 以下仅适用于离心轴流部分
      GPRC=GPRV/PKSI
      IF(GPRC.LT.GI(1)) GOTO 2
      IF(IIQP.GT.0) GOTO 4
      IQP=1
      IPE4=0
      QQS=PLA(4,M)
      PR=PLA(2,M)
      CALL HARIKA(GB,GKA,AN2,A,B,X,Y,SIG,ALF,NTB,IKBD,
     *IQP,IWRITE,AKGK,NNS,NNF,IIS,IIF,QQS,QQF,IPE4,KPAH,FKPAH,
     *DSIG,IZAP,LSR,LCP,CU,EP,PR,PB,BRANCH,OPTIMS,GKI,
     *POINTS,BASE,STORE,GAN,PLANE,I)
      IQP=0
      IPE4=1
      IIQP=IIQP+1
      QQS=0.
      PR=0.
      DO 5 IS=1,IZ
      CALL DECOD(STORE(167+8*(IS-1)),5,A4,4,B4,4)
      CALL DECOD(STORE(168+8*(IS-1)),6,A5,4,B5,4)
      CALL DECOD(STORE(171+8*(IS-1)),6,A6,4,B6,3)
      SRN(1,IS)=A4
      SRN(2,IS)=A6
      SRN(3,IS)=A5
    5 CONTINUE
    4 CONTINUE
      DO 14 JT=1,19
      IF(GPRC.GT.GI(JT).AND.GPRC.LE.GI(JT+1)) GOTO 11
      GOTO 14
   11 PIC(M)=((PII(JT+1)-PII(JT))*GPRC+PII(JT)*GI(JT+1)-PII(JT+1)*
     *GI(JT))/(GI(JT+1)-GI(JT))
      ETC(M)=((ETI(JT+1)-ETI(JT))*GPRC+ETI(JT)*GI(JT+1)-ETI(JT+1)*
     *GI(JT))/(GI(JT+1)-GI(JT))
   14 CONTINUE
      IF(ITC.EQ.1) GOTO 2
      IF(IDOP.NE.0) GOTO 12
      IF(M.EQ.1) GOTO 12
      PLAG=PLA(1,M-1)
      PLAP=PLA(2,M-1)
      PLAE=PLA(3,M-1)
      IDOP=0
   12 CONTINUE
      IF(ETC(M).NE.0..AND.IDOP.NE.0) ITC=1
      IF(ETC(M).NE.0.) GOTO 13
      PLA(1,M)=(PLA(1,M)+PLAG)/2.
      PLA(2,M)=(PLA(2,M)+PLAP)/2.
      PLA(3,M)=(PLA(3,M)+PLAE)/2.
      IDOP=IDOP+1
      GOTO 15         ! 重新计算？
   13 CONTINUE
      PIOC(I,M)=PLA(2,M)*PKSI*PIC(M)
      ETOC(I,M)=(PIOC(I,M)**.2857-1.)/((1.+(PLA(2,M)**.2857-1.)
     */PLA(3,M))*(1.+(PIC(M)**.2857-1.)/ETC(M))-1.)
      GPROC(I,M)=PLA(1,M)
      WRITE(7,113)GPROC(I,M),PIOC(I,M),ETOC(I,M),PLA(2,M),PLA(3,M) 
     *,PIC(M),ETC(M)
    2 CONTINUE
    3 CONTINUE
      WRITE(7,118)(SRN(1,IT),IT=1,IZ)  ! 输出
      WRITE(7,119)(SRN(2,IT),IT=1,IZ)
      WRITE(7,120)(SRN(3,IT),IT=1,IZ)
      WRITE(7,17)
      IF(IOCK.EQ.1) GOTO 7
      GOTO 10
    7 CONTINUE
      DO 9 M1=1,LB3
      GPROC(I,M1)=PLA(1,M1)
      PIOC(I,M1)=PLA(2,M1)
      ETOC(I,M1)=PLA(3,M1)
    9 CONTINUE
   10 CONTINUE
      DO 3559 J=1,IZ
      B(5,J)=B(5,J)-TETA(I,J+1)
      B(6,J)=B(6,J)-TETA(I,J+1)
 3559 CONTINUE
      GB(4)=GB(4)-TETA(I,1)
      GOTO 37
      ! 以上计算轴流部分
   36 CONTINUE
      ! 以下至37为离心压气机计算
      WRITE(7,114)AN2(1,1)
      WRITE(7,115)
      DO 38 J=1,20
      GPROC(I,J)=GI(J)
      PIOC(I,J)=PII(J)
      ETOC(I,J)=ETI(J)
      IF(GI(J).EQ.0.) GOTO 38
      WRITE(7,113)GI(J),PII(J),ETI(J)
   38 CONTINUE
      DO 39 K=1,IZ2
      TETA(I,K)=0.
   39 CONTINUE
      IZ=0
   37 CONTINUE
 3557 CONTINUE
      ! 输出结果
C      open(7010,FILE='bladevars.s12s1')
	open(7010,FIle='Par_BLADE.1DDtoANA')        ! 输出一维特性计算所需的叶片几何参数（用于向导模式估算）
      do l=1,iz
      WRITE(7010,'(20f15.6)')(A(J,L),J=1,20)  ! 动叶栅  
      enddo
      do l=1,iz
      WRITE(7010,'(20f15.6)')(B(J,L),J=1,20)  ! 静叶栅
      enddo
                     
      open(7011,file="Par_BC.1DDtoANA",status="old", position="append",   ! 输出一维特性计算所需边界条件
     *action="write")
      write(7011,*)A(1,1)    ! 进口外径
      close(7011)
      close(7011)
      close(7010)
      close(7002)
      close(7006)
      IF(IOCK.EQ.1.AND.IR.NE.2) PIK=PK
      IF(IOCK.EQ.2) GPRR=GR
      CALL GH(NV,GPROC,PIOC,ETOC,TETA,IZ,AN,PIK,GPRR,IOCK)
   33 CONTINUE
   17 FORMAT(80(1H-))
  113 FORMAT(7F11.4)
  114 FORMAT(3X,'characteristic centrifugal stage N OTN=',F5.3)
  115 FORMAT(5X,2HGC,10X,3HPIC,10X,3HETC)
  116 FORMAT(5X,4HG PR,6X,4HPI K,9X,3HETA,7X,4HPI A,7X,5HETA A,6X,
     *4HPI C,7X,5HETA C)
  117 FORMAT(3X,'total characteristic',5X,'N =',F5.3)
  118 FORMAT(3X,'NU TEK=',10F8.4)
  119 FORMAT(3X,'NU S K=',10F8.4)
  120 FORMAT(3X,'NU S S=',10F8.4)
      stop
      END
      SUBROUTINE BBOD1(OB,PK,G3,PBX,T1,XTA3,AKP,IZ,K,KF,KDB,UBX,PBHA,
     *D3,D31,ID1,IZ2,BXD1,ALF1,IOCK,N,AN,TETA,GB,PSI)
      DIMENSION D3(9),D31(30),ID1(33),BXD1(10),
     *C1AS(30),C2AS(30),AN(10,3),TETA(10,30),GB(10),YK(30)
      REAL ID1
      COMMON/YKC/YK
      COMMON/AH11/AH1(30),CMKA,CATAK
      COMMON/C12AS/C1AS,C2AS
      COMMON/BIS/BISR,BISS                             ! 级间引气参数
      DIMENSION BISR(30),BISS(30)
      DIMENSION AIB(100)
      READ(4,*)PBX,T1,UBX,PBHA,ID1(4),ID1(32)
      write(0,*)PBX,T1,UBX,PBHA,ID1(4),ID1(32)
      READ(4,*)OB,PK,G3,XTA3,ID1(13),AKP
      write(0,*)OB,PK,G3,XTA3,ID1(13),AKP
      READ(4,*)IZ,KF,IOCK,CMKA,CATAK
      write(0,*)IZ,KF,IOCK,CMKA,CATAK
      K=1
      KDB=1
      WRITE(7,2)
      WRITE(7,17)
      WRITE(7,5)
      open(7011,file="Par_BC.1DDtoANA")
      
      WRITE(7,6)PBX,T1,UBX,PBHA,ID1(4),ID1(32)
      write(7011,*)pbx    ! 进口总压
      write(7011,*)T1     ! 进口总温
      write(7011,*)id1(4) ! 进口气流角

     
      write(7006,*)'TPG'               ! 标识1
      write(7006,'(2f15.8)') PBX,T1    ! 输出进口总压和总温
      IZ2=IZ+1                            ! 用于标识每级间隙的个数，即IZ+1
      IZ3=IZ*2+1                          ! 用于标识每个叶栅的个数，即IZ*2+1
      WRITE(7,17)
      WRITE(7,10)
      WRITE(7,57)OB,PK,G3,XTA3,ID1(13),AKP
      write(7011,*)ob  ! 转速
      write(7011,*)pk  ! 压比
      write(7011,*)g3  ! 流量
      
      write(7006,'(2f15.8,i)') G3,OB,iz ! 输出设计流量，转速，级数
      WRITE(7,17)
      DO 27 I=1,9
      D3(I)=0.            ! 初始化D3数组为0
   27 CONTINUE
      WRITE(7,12)
      WRITE(7,26)IZ,KF,IOCK,CMKA,CATAK
      write(7011,*)iz   ! 级数
      write(7011,*)kf   ! 叶型编码
      write(7011,*)iock ! 压气机形式编码
      
      close(7011)
      READ(4,*)(D31(I),I=1,IZ2)   ! 读取每级入口处工作轮外径和压气机出口外径
      WRITE(7,17)
      ! TODO: 此处if语句可优化
      IF(D31(IZ2).EQ.0.) GOTO 28  ! 当压气机出口外径为0时，数据错误，终止程序
      GOTO 29
   28 PRINT*,'PUT DK1(IZ+1)'
      STOP
   29 CONTINUE
      J=IZ2/2                 ! J为？
      IF(D31(2).NE.0) GOTO 40
      IF(IZ.EQ.1) GOTO 40
      DO 41 I=1,J             ! 当级数不为1且没有给出其他级外径时，自动计算每级外径
   41 D31(I)=D31(1)
      J1=IZ-J
      DO 42 I=1,J1
   42 D31(I+J)=D31(1)-(D31(1)-D31(IZ2))/(IZ2-J)*I
   40 CONTINUE
      IF(IZ2.GT.10)GOTO 33    ! 每行输出10个级的外径
      WRITE(7,23)(D31(I),I=1,IZ2)
      WRITE(7,17)
      GOTO 34
   33 WRITE(7,23)(D31(I),I=1,10)
      WRITE(7,17)
      WRITE(7,23)(D31(I),I=11,IZ2)
      WRITE(7,17)
      IF(IZ2.GT.20)GOTO 336
      GOTO 34
  336 WRITE(7,23)(D31(I),I=21,IZ2)
      WRITE(7,17)
   34 CONTINUE
      READ(4,*)(C1AS(I),I=1,IZ2)  ! 读取工作轮前轴向速度
      write(0,*)(C1AS(I),I=1,IZ2)
      WRITE(7,35)(C1AS(I),I=1,IZ2)
      WRITE(7,17)
      READ(4,*)(C2AS(I),I=1,IZ)   ! 读取工作轮后轴向速度
      write(0,*)(C2AS(I),I=1,IZ)
      WRITE(7,36)(C2AS(I),I=1,IZ)
      WRITE(7,17)
      ID1(1)=C1AS(1)
      ID1(2)=(C1AS(1)+C1AS(IZ2))/2.
      ID1(3)=C1AS(IZ2)
      READ(4,*)(AH1(I),I=1,IZ)    ! 读取每一级压头损耗系数
      write(0,*)(AH1(I),I=1,IZ)
      WRITE(7,37)(AH1(I),I=1,IZ)
      WRITE(7,17)
      ID1(7)=AH1(1)
      ID1(9)=AH1(IZ)
      ID1(8)=AH1(1)
      READ(4,*)(YK(I),I=1,IZ)     ! 读取叶冠的展弦比变化系数
      write(0,*)(YK(I),I=1,IZ)
      WRITE(7,38)(YK(I),I=1,IZ)
      WRITE(7,17)
      READ(4,*)ID1(5),ID1(6),ID1(10),ID1(11),ID1(12)              ! 读取反力度和理论压头减小系数相关数据
      write(0,*)ID1(5),ID1(6),ID1(10),ID1(11),ID1(12)
      READ(4,*)ID1(14),ID1(15),ID1(16),ID1(17),ID1(18),ID1(24)    ! 读取工作轮和导向器叶片的相关数据
      write(0,*)ID1(14),ID1(15),ID1(16),ID1(17),ID1(18),ID1(24)
      READ(4,*)ID1(19),ID1(20),ID1(21),ID1(27),ID1(25),ID1(26)    ! 读取入口导叶的相关数据
      write(0,*)ID1(19),ID1(20),ID1(21),ID1(27),ID1(25),ID1(26)
      READ(4,*)BXD1(1),BXD1(2),BXD1(3),ALF1,BXD1(8),PSI           ! 读取离心级的相关数据
      IF(PSI.LT.0..OR.PSI.GT.60.) PSI=0.                          ! 工作轮出口通道中线与转轴法线的夹角应处于0~60度范围内
      write(0,*)BXD1(1),BXD1(2),BXD1(3),ALF1,BXD1(8),PSI
      READ(4,*)BXD1(6),BXD1(7),BXD1(4),BXD1(9),BXD1(10)           ! 读取离心级的相关数据
      write(0,*)BXD1(6),BXD1(7),BXD1(4),BXD1(9),BXD1(10)
      READ(4,*)GB(8),GB(9),GB(10),N                               ! 读取流量特性线计算指标
      write(0,*)GB(8),GB(9),GB(10),N
      DO 79 I=1,N
      READ(4,*)(AN(I,J),J=1,2)                                    ! 读取特性线的转速百分比及流量系数计算初始点
      write(0,*)(AN(I,J),J=1,2)
   79 CONTINUE
      DO 77 I=1,N
      READ(4,*)(TETA(I,J),J=1,IZ2)                                ! 读取可调叶冠转动角数组
      write(0,*)(TETA(I,J),J=1,IZ2)
   77 CONTINUE
      ! 添加级间引气参数
      DO 814 I=1,IZ
      READ(4,*)BISR(I),BISS(I)
  814 CONTINUE
      ID1(22)=0.                                                  ! 未知用途
      ID1(23)=0.
      ID1(28)=4300.
      ID1(29)=4300.
      ID1(30)=4300.
      ID1(31)=100.
      ID1(33)=1.
      WRITE(7,8)
      WRITE(7,3)ID1(5),ID1(6),ID1(10),ID1(11),ID1(12)
      WRITE(7,17)
      WRITE(7,19)
      WRITE(7,22)ID1(14),ID1(15),ID1(16),ID1(17),ID1(18),ID1(24)
      WRITE(7,17)
      WRITE(7,24)
      WRITE(7,21)ID1(19),ID1(20),ID1(21),ID1(27),ID1(25),ID1(26)
      WRITE(7,17)
      WRITE(7,13)
      WRITE(7,16)BXD1(1),BXD1(2),BXD1(3),ALF1,BXD1(8),PSI
      WRITE(7,17)
      WRITE(7,20)
      WRITE(7,16)BXD1(6),BXD1(7),BXD1(4),BXD1(9),BXD1(10)
      WRITE(7,17)
      WRITE(7,76)
      WRITE(7,75)GB(8),GB(9),GB(10),N
      WRITE(7,17)
      WRITE(7,74)
      DO 70 I=1,N
      WRITE(7,73)(AN(I,J),J=1,2)
   70 CONTINUE
      WRITE(7,17)
      WRITE(7,72)
      DO 71 I=1,N
      WRITE(7,78)(TETA(I,J),J=1,IZ2)
   71 CONTINUE
      WRITE(7,17)
   74 FORMAT(6X,1HN,7X,1HQ)
   73 FORMAT(10F8.4)
   78 FORMAT(10F8.2)
   72 FORMAT(6X,4HTETA)
   75 FORMAT(2X,F8.6,2X,F8.6,2X,F8.6,5X,I2)
   76 FORMAT(6X,2HDQ,8X,2HEQ,8X,2HES,7X,2HNV)
    1 FORMAT(8F9.1)
    4 FORMAT(6I9)
    2 FORMAT(35X,'计算结果输出列表')
    3 FORMAT(5(3X,F5.3))
   26 FORMAT (5X,I2,10X,I1,10X,I2,8X,F4.2,8X,F4.2)
   57 FORMAT(2X,F8.0,5(3X,F7.3))
    5 FORMAT (5X,2HP1,5X,2HT1,6X,2HSP,6X,4HSVNA,4X,4HALF1,3X,4HALFZ)
    6 FORMAT (F8.1,2X,F6.1,2X,F6.3,
     *2X,F6.3,2X,F6.2,2X,F6.2)
    8 FORMAT (5X,4HTAUS,3X,4HDTAU,5X,3HKH1,5X,3HDKH,5X,5HKHMIN)
   10 FORMAT (5X,1HN,9X,2HPK,9X,1HG,9X,3HETA,9X,2HKG,6X,3HKPI)
   12 FORMAT (8X,1HZ,8X,2HKF,8X,4HIOCK,7X,4HCMKA,7X,5HCATAK)
   20 FORMAT (7X,3HD1P,6X,3HD1V,7X,3HBL2,6X,4HSKOL,6X,4HOTBR)
   13 FORMAT (7X,3HPIS,6X,3HZST,7X,3HALP,6X,4HALF1,6X,4HAL1C,6X,3HPSI)
   16 FORMAT (6(4X,F6.3))
   17 FORMAT (80(1H-))
   19 FORMAT (5X,4HHRK1,5X,4HHRKZ,5X,5HXFVNA,4X,4HXFRK,5X,4HXFNA,
     *5X,4HPRK1)
   35 FORMAT (6X,3HC1A/10F8.1)
   36 FORMAT (6X,3HC2A/10F8.1)
   37 FORMAT (6X,2HHZ/10F8.4)
   38 FORMAT (6X,2HYK/10F8.2)
   23 FORMAT (6X,3HDK1/10F8.4)
   22 FORMAT (6(4X,F5.3))
   24 FORMAT(5X,4HPVNA,5X,4HHVNA,4X,6HB/TVNA,4X,5HCMVNA,4X,3HLOS,4X,
     *4HDLOS)
   21 FORMAT (6(4X,F5.3))
      RETURN
      END
      SUBROUTINE BBOD2(IOCK,PIK,GPRR,GB,GKA,AN,A,B,IZ,KF,NV,TETA,OBC,
     *ETCBK,D2C,BET2L,ALA3)
      COMMON/BIS/BIS(30)
      DIMENSION GB(10),GKA(4,30),AN(10,3),A(20,30),B(20,30),TETA(10,30)
      READ(4,*)IOCK,PIK,GPRR
      WRITE(7,4)
      IF(IOCK.EQ.1) WRITE(7,1)
    1 FORMAT(5X,'output axi.')
      IF(IOCK.EQ.2) WRITE(7,2)
    2 FORMAT(5X,'output cen.')
      WRITE(7,4)
      WRITE(7,50)PIK,GPRR
   50 FORMAT(5X,'PIR=',F6.3,5X,'GPRR=',F6.2)
      IF(IOCK.EQ.2) GOTO 3
      WRITE(7,4)
    4 FORMAT(80(1H-))
      READ(4,*)IZ,NV,KF
      WRITE(7,6)IZ,NV,KF
    6 FORMAT(5X,'IZ=',I2,5X,'NV=',I2,5X,'KF=',I2)
      WRITE(7,4)
      READ(4,*)(GB(I),I=1,10)  
      WRITE(7,7)
    7 FORMAT(5X,2HP1,5X,2HT1,4X,3HUK1,4X,4HALF1,10X,1HR,6X,1HK,6X,2HDQ,
     *6X,2HEQ,6X,2HES)
      WRITE(7,8)(GB(I),I=1,10)
    8 FORMAT(F8.1,1X,F6.1,2X,F6.2,1X,F5.1,F5.2,F8.2,F8.4,F8.6,F8.6,F8.6)
      WRITE(7,4)
      DO 9 J=1,IZ
      READ(4,*)(A(I,J),I=1,20)
    9 CONTINUE
      WRITE(7,12)
      DO 10 J=1,IZ
      WRITE(7,11)(A(I,J),I=1,20)
   10 CONTINUE
   11 FORMAT(10F8.4)
   12 FORMAT(6X,'寑憫垈 A(I,J)')
      WRITE(7,4)
      DO 13 J=1,IZ
      READ(4,*)(B(I,J),I=1,20)
   13 CONTINUE
      WRITE(7,16)
      DO 14 J=1,IZ
      WRITE(7,15)(B(I,J),I=1,20)
   14 CONTINUE
   15 FORMAT(10F8.4)
   16 FORMAT(6X,'寑憫垈 B(I,J)')
      WRITE(7,4)
      DO 60 J=1,IZ
      READ(4,*)(GKA(I,J),I=1,4)
   60 CONTINUE
      WRITE(7,61)
      DO 62 J=1,IZ
      WRITE(7,63)(GKA(I,J),I=1,4)
   62 CONTINUE
   61 FORMAT(6X,'GKA')
   63 FORMAT(4F8.4)
      WRITE(7,4)
      DO 79 I=1,NV
      READ(4,*)(AN(I,J),J=1,2)
   79 CONTINUE
      WRITE(7,74)
      DO 70 I=1,NV
      WRITE(7,73)(AN(I,J),J=1,2)
   70 CONTINUE
   74 FORMAT(6X,1HN,7X,1HQ)
   73 FORMAT(10F8.4)
      WRITE(7,4)
      IZ2=IZ+1
      DO 77 I=1,NV
      READ(4,*)(TETA(I,J),J=1,IZ2)
   77 CONTINUE
      WRITE(7,72)
      DO 71 I=1,NV
      WRITE(7,78)(TETA(I,J),J=1,IZ2)
   71 CONTINUE
   72 FORMAT(6X,4HTETA)
   78 FORMAT(10F8.2)
      WRITE(7,4)
      ! 添加级间引气参数
      READ(4,*)(BIS(I),I=1,IZ)
      GOTO 40
    3 CONTINUE
      READ(4,*)OBC,ETCBK,NV
      WRITE(7,30)OBC,ETCBK,NV
   30 FORMAT(5X,'NPR=',F8.1,5X,'ETR=',F5.3,5X,'NV=',I2)
      WRITE(7,4)
      READ(4,*)D2C,BET2L,ALA3
      WRITE(7,31)D2C,BET2L,ALA3
   31 FORMAT(6X,'D2=',F6.4,4X,'BET2L=',F6.2,4X,'LAM3=',F6.4)
      WRITE(7,4)
      READ(4,*)(AN(I,1),I=1,NV)
      WRITE(7,32)
   32 FORMAT(5X,7HAN(I,1))
      WRITE(7,73)(AN(I,1),I=1,NV)
      WRITE(7,4)
   40 CONTINUE
      RETURN
      END
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
      ! ---   IKBD    压气机级序号，从此级开始第二级???
      COMMON/A1/D1/A2/R1/A3/YY,H1,H2,T1,T2,CK
      ! A1
      ! ---   D1      转子进口处轮毂比d1
      ! A2
      ! ---   R1      ???
      ! A3
      ! ---   YY      工作轮叶片弦长
      ! ---   H1      ???
      ! ---   H2      ???
      ! ---   T1      动叶进口节距
      ! ---   T2      ???
      ! ---   CK      ???
      COMMON/A5/HRK
      ! A5
      ! ---   HRK     ???
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
      ! ---   UK1     转子进口周向速度
      ! A8
      ! ---   ALFA1   工作轮进口处气流角
      ! A9
      ! ---   RRR     气体常数R
      ! ---   AT1     转子进口总温
      ! ---   UK2     转子出口周向速度
      ! ---   UK4     静子出口周向速度
      ! A10
      ! ---   AKH     ???
      COMMON/A11/C2AO,C4AO,ALFA4O/A12/AKP/A13/AK
      ! A11
      ! ---   C2AO    ???
      ! ---   C2AO    ???
      ! ---   ALFA4O  导向器出口处最佳气流角
      ! A12
      ! ---   AKP     ???
      ! A13
      ! ---   AK      ???
      COMMON/A14/HB3,HB4,CA,BT3,BT4
      ! A14
      ! ---   HB3     ???
      ! ---   HB4     ???
      ! ---   AT1     ???
      ! ---   BT3     ???
      ! ---   BT4     ???
      COMMON/A15/AKR
      ! A15
      ! ---   AKR     ???
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
      ! ---   PPI     ???
      COMMON/A29/AT2
      ! A29
      ! ---   AT2     ???
      COMMON/A30/C1A,AML1
      ! A30
      ! ---   C1A     ???
      ! ---   AML1    ???
      COMMON/A31/F1,XFK,XFA/A32/ANT,AKG
      ! A31
      ! ---   F1      ???
      ! ---   XFK     ???
      ! ---   XFA     ???
      ! A32
      ! ---   ANT     ???
      ! ---   AKG     ???
      COMMON/A33/QLA1C,HQCK,HQC
      ! A33
      ! ---   QLA1C   ???
      ! ---   HQCK    ???
      ! ---   HQC     ???
      COMMON/A34/PM
      ! A34
      ! ---   PM      第一级入口总压???
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
     *6X,4HPI K,9X,3HKPD,9X,1HG,10X,4HG PP,6X,6HGPRBIX)
  117 FORMAT(5X,'characteristic axi. stage',5X,'N =',F5.3,
     *5X,'UK =',F5.1)
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
      GB(1)=GB(1)/9.81        ! 第一级入口总压
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
      DO L=IIS,IIF
          IF(GKA(2,L).EQ.0.) GKA(2,L)=1.
          IF(GKA(3,L).EQ.0.) GKA(3,L)=0.98
          IF(GKA(4,L).EQ.0.) GKA(4,L)=1.             !流量、压比、效率的修正系数，及参考值；落后角修正系数在下面
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
      DO L=1,81               ! 这个循环是啥意思
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
      ! /************ 级循环，行号 1 ************/
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
      IF(UKM.LE.150.AND.RK(16).NE.0.) THEN
          RK(10)=0   !!!!!!!!!!!!!!
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
      IF(I.LE.IIS) THEN
          ! 第一级压气机的参数处理
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
      IF(FORM2.GT.0.000001) THEN
          GB(1)=GB(1)*FORM2
          STAGE(22)=STAGE(22)*FORM2
      ENDIF
      IF(INDEX4.LT.1.AND.I.LE.IIS) THEN
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
      ! /************ 级循环结束，行号 1 ************/
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
      GPRBIX=G*10328.746/GB(1)*SQRT(GB(2)/288.)       !??????
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
              WRITE(7,117)BN,CN
              write(numberchar,'(f5.3)') bn
              if(outflag.eq.1) then
              write(7002,*)'zone t="character'//trim(numberchar)//'"'
              endif
          endif
          IF(KRAN(39).NE.0.OR.KRAN(40).NE.0) THEN
              WRITE(7,116)
              WRITE(7,113)QL1,PP,AKPX,G,GPP,GPRBIX    !输出
          ELSE
              IF(NG-1.EQ.((NG-1)/10)*10) WRITE(7,116)
              WRITE(7,113)QL1,PP,AKPX,G,GPP,GPRBIX    ! 输出流量系数、等熵效率等参数
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
      GAN(2)=GPP      ! 流量
      GAN(3)=QL1      ! 流量系数
      GAN(4)=PP       ! PP???
      GAN(5)=AKPX     ! 等熵效率
      GAU=GB(1)*9.81  ! 总压x9.81
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
  198 FORMAT('    媹崍?悁亷棃?悈唸寧?媴倕?儛€崍枦 搼拵墬垈帒拡')
      RETURN
      END
      SUBROUTINE SFROMT(T,S)
      T1=T/100.
      S=SQRT((0.210544*T1+3.55644)*T1-2.05832)
      S=1000.*(-0.4114208*T1+4.70161+S)
      RETURN
      END
      SUBROUTINE TFROMS(S,T)
      S1=S/1000.
      T=SQRT((1.23576*S1-20.2081)*S1+86.7518)
      T=1000.*(0.996738*S1-8.99432+T)
      RETURN
      END
      SUBROUTINE IFROMT(T,AI)
      T1=T/100.
      AI=SQRT((96.58115*T1-1393.581)*T1+6252.329)
      AI=1000.*(108.0157*T1-77.5678+AI)
      RETURN
      END
      SUBROUTINE AKPDI(T1,T2,PI,R,AKPD)
      CALL IFROMT(T1,AI1)
      CALL SFROMT(T1,S1)
      SAD=S1+R*ALOG(PI)
      CALL TFROMS(SAD,TAD)
      CALL IFROMT(T2,AI2)
      CALL IFROMT(TAD,AIAD)
      AKPD=(AIAD-AI1)/(AI2-AI1)
      RETURN
      END
      SUBROUTINE DLIHA(B,BB,D1,TGT,S,SB,DL,DB,K,IZ,DS)
      DIMENSION B(30,2),BB(30,2),D1(30,2),TGT(30,2),SI(30),
     *DI(30,2)
      COMMON/DL1/DLBHA,DLPEP,SI,DI
      DBT=DL
      DL=0.
      DO 4 I=1,IZ
      DO 4 J=1,2
      AJ=I
        AJ1=AJ-1.
        Z=IZ
        AZ=Z-1.
      IF(IZ.EQ.1)SI(1)=S
      IF(IZ.GE.2)SI(I)=S+DS*AJ1**.5/AZ**.5  !轴向间隙
      SAS=SI(I)*B(I,1)
      IF(SAS.LT..005)SAS=.005
      SI(I)=SAS/B(I,1)
      TET=TGT(I,J)
      P=1-D1(I,J)
      IF(J.EQ.2)GOTO 1
      O=(1.+.3*P)/(1.+.75*P*(BB(I,J)-1.))
        GOTO 2
    1 O=(1.+.2*P)/(1.+.5*P*(BB(I,J)-1.))
    2 IF(J.EQ.2.AND.I.EQ.IZ)GOTO 5
      DL1=B(I,J)*(SIN(TET)*O+SI(I)) !转子轴向长度
      GOTO 3
    5 DL1=B(I,J)*SIN(TET)*O
    3 DL=DL+DL1
      DI(I,J)=DL1
      SI(I)=SI(I)*B(I,1)
    4 CONTINUE
      IF(DB.LT..001)GOTO 9
      DBT=DBT/DB
      IF(DBT.GT.1.)DBT=1.
      IF(K.EQ.6.OR.K.EQ.7.OR.K.EQ.3)
     *DP=DB*(.585*(1.-DBT)+.35)
      IF(K.LT.6.AND.K.NE.3)DP=0.
      IF(K.GT.7)
     *DP=DB*(.725*(1.-DBT)+.185)
      GOTO 17
    9 DP=0.
   17 DLBHA=DL+SB
        DLPEP=DP
      DL=DL+DP+SB
      IF(K.EQ.1)DL=DLBHA
      SB=DP
    8 FORMAT (5F10.4)
      RETURN
        END
      SUBROUTINE BEC(DP,DBT,R0,BKA,BT,CKA,UK1,IZ,KDB,ID1,
     *PCYM,GBHA,PEP,BBHA,T)
      COMMON/BE/GCYM,GKOMP,GPEP
      DIMENSIONDP(30,2),DBT(30,2),R0(30,2),BKA(30,2),BT(30,2),
     *CKA(30,2),ID1(33),T(30)
      REAL ID1
      DATAP/3.14159/,PA/57.295828/
      CYM1=0.
        CYM2=0.
        CYM3=0.
        CYM4=0.
      CT=7850.
      CIGP=ID1(31)
      PP=ID1(28)
        PC=ID1(29)
        PL=ID1(30)
      DHPEP=ID1(23)
      PCYM=PCYM*ID1(33)
      IZ7=IZ+1
      GPEP=410.*PEP*DHPEP**2
      GOTO (2,2,2,4,6,6,6,8 ),KDB
    2 A=.00088*UK1+.17*DBT(1,1)-.252
      B=.003109*UK1/100.-.000864
      GOTO 10
   14 A=1.52-2.6*DBT(IZ7,1)*DP(IZ7,1)/DP(1,1)
    7 B=.00833*UK1/100.-.02483
      GOTO 10
    6 A=.555-.0028*CIGP
      B=.000454*PCYM+.020
      GOTO 10
    8  A=.555-.0028*CIGP
      B=.001375*PCYM-.0161
   10 DO 15 I=1,IZ
      IF(T(I  )-630.)17,17,19
   19 PC=CT
      IF(T(I)-680.)17,17,20
   20 PA=CT
      IF(T(I)-790.)17,17,18
   18 PL=CT
        PP=CT
   17 B=BKA(I,1)
        BA=BKA(I,2)
      I4=I+1
      DK1=DP(I,1)
        DK2=DP(I,2)
        DK4=DP(I4,1)
      DB1=DBT(I,1)*DK1
      DB2=DBT(I,2)*DK2
      DB4=DBT(I4,1)*DK4
      BTK=BT(I,1)
        BTA=BT(I,2)
      R1=R0(I,1)
        R2=R0(I,2)
        R4=R0(I4,1)
      CK=CKA(I,1)
        CA=CKA(I,2)
      H1=(DK1-DB1)/2.
        H2=(DK2-DB2)/2.
        H4=(DK4-DB4)/2.
        DPK=(DK1*R1+DK2*R2)/2.
      DHA=(DK2*R2+DK4*R4)/2.
      HPK=(H1+H2)/2.
        HHA=(H2+H4)/2.
      BB=B+BA
      IF(I.NE.1)GOTO 12
      IF(GBHA.GT..001)GOTO 12
      BB=BB+BBHA
        BA=BA+BBHA
   12 CYM1=PP*BB*(DB1**2+DB4**2+DB1*DB4)+CYM1
      CYM2=PC*BB*(DK1+DK4)+CYM2
      CYM3=Pl*DPK*BTK*HPK*B*CK+CYM3
      CYM4=PA*DHA*BTA*HHA*BA*CA+CYM4
   15 CONTINUE
      GKOMP=P*((A*CYM1)/12.+B*(DP(1,1)+DP(IZ7,1))/4*CYM2+
     *CYM3+CYM4)+GBHA
      GCYM=GKOMP+GPEP
      GOTO 27
    4 IF(DBT(IZ7,1).LT..55)GOTO 14
      A=1.52-1.43*DP(IZ7,1)/DP(1,1)
        GOTO 7
   27 RETURN
        END
      REAL FUNCTION XHA(X1)               !由温度求焓
        X=X1/100.
      XHA=10.**3*(108.0157*X-77.5678+SQRT(96.58115*X**2
     *-1393.581*X+6252.329))
        RETURN
        END
      REAL FUNCTION DF(P,F,D)
        DF=SQRT(D**2+4.*F/P)
        RETURN
        END
      REAL FUNCTION F(CA,C,T,P,R,G,X)
      XC=X-C**2/2.                        !进口静焓
        TC=TXH(XC)                        !进口静温
      PC=P/EXP((XHP(T)-XHP(TC))/R)        !进口静压
      F=R*G*TC/PC/CA                      !进口面积
        RETURN
        END
      REAL FUNCTION TAN(X)
        TAN=SIN(X)/COS(X)
        RETURN
        END
      REAL FUNCTION CAPOF(F,T,P,R,G,X,A1)
      CA1=1.2*R*G*T/F/P
        SIA=SIN(A1)
    2 C1=CA1/SIA
        AI=X-C1**2/2.
        TC=TXH(AI)
      CA=R*G*TC/F/(P/EXP((XHP(T)-XHP(TC))/R))
      IF(ABS(CA-CA1).LE..1)GOTO 4
      CA1=CA
        GOTO 2
    4 CAPOF=CA
      RETURN
        END
      REAL FUNCTION FD(P,D,B)
        FD=P*(D**2-B**2)/4.
        RETURN
        END
      REAL FUNCTION RD(D)
        RD=SQRT((1.+D**2)/2.)
        RETURN
        END
      REAL FUNCTION TXP(X3)
        X=X3/1000.
      TXP=10.**3*(.996738*X-8.99432+SQRT(1.23576*X**2
     *-20.2081*X+86.7518))
        RETURN
        END
      REAL FUNCTION XHP(X2)             !熵函数
        X=X2/100.
      XHP=10.**3*(4.70161-.4114208*X+SQRT(.210544*
     *X**2+3.55644*X-2.05832))
        RETURN
        END
      REAL FUNCTION TXH(X3)
        X=X3/1000.
      TXH=.93351754*X+66.388971-SQRT(.0072137889*X**2
     *-10.124088*X+4611.07977)
        RETURN
        END
      REAL FUNCTION TKP(X4)
        X=X4/100.
      TKP=.41394-.0015178*X**4+.0371707*X**3+
     *.051916*X**2+82.5818*X
        RETURN
        END
      SUBROUTINE PCPM(OB1,PK1,G31,PBX1,T11,XTA31,AKP1,IZZ,
     *K1,KF1,KDB1,UBX1,PBHA1,AD3,AD31,AD32,AID1,AID2,
     *EA,EB,EGKA,EGB,KKP,ESI,EDI,EAA,CYMP1,IZ21)
      COMMON/MOD/IMODD(5),CA(15),HZ(15),ALF(15),RC2(15),AKPK(15),
     *AKPDC(15)
      DIMENSION AA(30,65),EAA(30,65)
      DIMENSION EA(20,30),EB(20,30),EGKA(4,30),EGB(10),CYMP(20),
     *CYMP1(20),AD3(9),AD31(30),AD32(30),AID1(33),AID2(24),KKP(30,3),
     *ESI(30),KP(30,3),SI(30),EDI(30,2),DI(30,2)
      DIMENSION A(20,30),B(20,30),GKA(4,30),GB(10)
      DIMENSION D3(9),D31(30), D32(30)
      REAL ID1(33),ID2(24)
      COMMON/BB/D3,D31,D32,ID1
      COMMON/BB2/IZ1,IZ2,IZ3
      COMMON/Z4/K,IV1
      COMMON/BEH/ID2
      COMMON/BB1/OB,PK,G3,PBX,T1,XTA3,AKP,IZ,  KF,KDB,UBX,PBHA
      COMMON/PR1/KP
      COMMON/Z19/GKA,GB
      COMMON/CB/A,B
      COMMON/CY/CYMP,AA
      COMMON/DL1/DLBHA,DLPEP,SI,DI
      IZ1=IZZ
      IZ2=IZ21
        IZ3=IZ21
      OB=OB1                                   !转子的物理转速
        PK=PK1                                 !单级增压比
        G3=G31                                 !物理空气流量
        PBX=PBX1                               !压气机入口壳体前的总压
        T1=T11                                 !压气机前空气总温
      XTA3=XTA31                               !等熵效率值，待求则XTA3=0
        AKP=AKP1                               !单级增压比(PK)的增大系数
        IZ=IZZ                                 !轴流压气机级数
        K=K1                                   !轴流级组的给定直径符号，BBOD1中定义为1
        KF=KF1            !亚音速级的初始叶片形状符号1--BC-10,2--C-4,3--NACA-65
      KDB=KDB1                                  !压气机图形符号，BBOD1中定义为1
        UBX=UBX1                                !入口壳体的总压恢复系数
        PBHA=PBHA1                            !入口导向器的总压恢复系数
      DO 2 I=1,9
    2 D3(I)=AD3(I)                              !不用数组
    4 IF(K.EQ.0)GOTO 5
      IZ7=IZ+1
      DO 3 I=1,IZ7
      D31(I)=AD31(I)       !每级入口处工作轮盘的外径和压气机出口上的外径
    1 IF(K.LE.3)GOTO 3
      D32(I)=AD32(I)
    3 CONTINUE
    5 DO 7 I=1,33
    7 ID1(I)=AID1(I)
   16 IF(KDB.LT.9.)GOTO 9
      DO 8 I=1,24
    8 ID2(I)=AID2(I)
    9 CONTINUE
      CALL PCP1           ! 这是啥
      DO 20 I=13,16
      DO 20 J=1,IZ1
      A(I,J)=0.
   20 B(I,J)=0.
      DO 21 I=18,20
      DO 21J=1,IZ1
      A(I,J)=0.
       B(I,J)=0.
      B(8,J)=0.
   21 B(10,J)=0.
      GB(7)=0.
      GB(5)=0.
      GB(8)=EGB(8)
      GB(9)=EGB(9)
      GB(10)=EGB(10)
      DO 535 M=1,30
      GKA(1,M)=0.
      GKA(4,M)=0.
535   CONTINUE
      DO 10 I=1,10
   10 EGB(I)=GB(I)
      DO 11 I=1,20
      DO 11 J=1,IZ
      EA(I,J)=A(I,J)
        EB(I,J)=B(I,J)
   11 CONTINUE
      DO 12 I=1,4
      DO 12 J=1,IZ
   12 EGKA(I,J)=GKA(I,J)
      DO 14 I=1,IZ
      ESI(I)=SI(I)
      DO 14 J=1,2
      EDI(I,J)=DI(I,J)
   14 CONTINUE
      DO 19 I=1,20
   19 CYMP1(I)=CYMP(I)
      DO 15 I=1,30
      DO 15 J=1,65
   15 EAA(I,J)=AA(I,J)
      DO 17 I=1,IZ
      DO 17 J=1,3
   17 KKP(I,J)=KP(I,J)
      IZZ=IZ
      UBX1=UBX
      RETURN
        END
      SUBROUTINE PCP1
      DIMENSION PG(30,2),PW(30,2),CIG(9)
      DIMENSION D3(9),D31(30), D32(30)
      REAL ID1(33)
      DIMENSION SI(30),DI(30,2),CYMP(20),AA(30,65)
      DIMENSION CU(30,2),TE(30),AL1(30),UK(30)
      DIMENSION A(20,30),B(20,30),GKA(4,30),GB(10)
      DIMENSION DP(30,2),DC(30,2),DB(30,2),DBT(30,2),
     *R0(30,2),F0(30,2),CAI(30,2),BT(30,2),BE1(30),BE2(30),
     *AL2(30),AL4(30),DL(30,2),L(30,2),HZI(30),HAD(30),
     *PC(30,2),XT(30,2),ALX(30),YD(30,2),IZKA(30,2),HTI(30),
     *BKA(30,2),XPC(30,2),AI(30,2),DELT(30,2),Q(30,2),AGA1(30,2)
      REAL L,L1,L2,LBYX
      REAL*8 PI
      EXTERNAL XHA,XHP,TXP
      EXTERNAL TXH,TKP
      COMMON/AH11/AH1(30),CMKA,CATAK
      COMMON/Z8/DP,DC,DB
      COMMON/Z10/R0,F0,DBT
      COMMON/BB/D3,D31,D32,ID1
      COMMON/BB1/OB,PK,G3,PBX,T1,XTA3,AKP,IZ,  KF,KDB,UBX,PBHA
      COMMON/BB2/IZ1,IZ2,IZ3
      COMMON/Z19/GKA,GB
      COMMON/CYM/PCYM,XTCYM,DLKOMP,LBYX,HB,HB1,X
      COMMON/CY/CYMP,AA
      COMMON/BE/GCYM,GKOMP,GPEP
      COMMON/DL1/DLBHA,DLPEP,SI,DI
      COMMON/Z1/DP1,DC1,DB1,DPB,DCB,DBB/Z2/HZ1,HZC,HZB,HZK,
     *UK1,UKB/Z3/C1A,CAC,CAB/CB/A,B/Z4/K,IV1
      COMMON/Z5/GP,PA,P,R/Z13/BC3
      COMMON/Z7/TAY,TAYC,TAYB
      COMMON/Z17/FB
      COMMON/Z11/BE1,BE2
      COMMON/PP/PG,PW,CIG
      COMMON/Z12/BT
      COMMON/BX/ALF1,C1U,P1,CA,R1
      COMMON/Z15/YD,AlX/Z16/IZKA,BKA
      COMMON/PK/CAI,AL2,AL4,AL1,DL,L,HZI,HAD,CU,TE,UK,PC,XT,Q,AI
      COMMON/Z21/K10,K11,IDZ
      COMMON/PK1/ALFB
      COMMON/Z42/ABT,KC3
      COMMON/MOD/IMODD(5),CAK(15),HZE(15),ALFK(15),RC2(15),AKPK(15),
     *AKPDS(15)
      BHAL=0
      BHAZ=0
      KBX=0
      K10=1
      K11=0
      IDZ=0
      TAY=0.
      XTA=.85
      AKGB=1.
      P=3.14159
        PA=180./P
        C1A=ID1(1)                 !工作轮入口处气流的轴向速度
        IF(IMODD(2).EQ.1)C1A=CAK(1)
        R=287.023
              AKG=ID1(13)         !空气流量给定值的增大系数
        GP=G3*AKG                 !实际给定的空气流量
      BC3=ID1(24)                 !第一超音速级工作轮叶片的受风度
        DBMAX=ID1(22)             !取值为0
      SKB=ID1(25)                 !第一级相对轴向间隙
      DS=ID1(26)                  !最后一级和第一级相对轴向间隙
      HK1=ID1(10)                 !第一级理论压头减小系数
        DKH=ID1(11)               !依次每级中理论压头减小系数减小的步长
        HKMI=ID1(12)              !理论压头减小系数的最小值
      ALF1=ID1(4)/PA              !第一级工作轮盘前绝对气流的进气角
      TAYC=ID1(5)                 !在中径上，轴流第N级的反动度
        TAYB=ID1(6)               !在每一级中i>N，反动度变化的步长
      HVB=ID1(9)                  !末级压头损失系数
      IF(IMODD(1).EQ.1)HVB=HZE(IZ)
      HZC=ID1(8)                  !第一级压头损失系数
        CAC=ID1(2)   !第一级工作轮前进口轴向速度与压气机出轴向速度平均值
      CAB=ID1(3)                  !压气机出口轴向分速度
      ALFB=ID1(32)/PA             !末级动叶出口气流角，弧度制
      XA=XHA(T1)                  !有温度求焓
   11 IF(KDB.GE.6.AND.KDB.LT.9)GOTO 13
      IF(KDB.EQ.3)GOTO 13
      GOTO 12
   13 IF(UBX.GT.0.)GOTO 12
      KBX=1
        UBX=0.99
   12 P1=PBX*PBHA*UBX                   !第一级转子进口前的总压
      CALL PCPBX                        !轴流级进口计算
   15 IF(KBX.EQ.0)GOTO14
      KBX=0
      UBX=.025*DB1/DBMAX+.97
      GOTO 12
   14 IF(KDB.GT.8)GOTO 120
      UK(1)=UK1
      CU(1,1)=C1U
      AL1(1)=ALF1*PA
   50 PKP=PK*AKP
   60 PB=PKP*PBX
      IF(XTA3.GT..01)XTA=XTA3
      DO 108 I=1,IZ
      UK(I)=UK1*D31(I)/D31(1)
      A(20,I)=AH1(I)*UK(I)**2
  108 CONTINUE
   63 CALL PCPBYX(PKP,PB,ALFB,AKGB,DBTB,TB,XA,XTA)   !轴流级出口计算
        UKB=UK1*DPB/DP1
      HZB=HVB*UKB**2
      IF(IZ.EQ.1)HZ1=ID1(7)*UK1**2
      CALL ZC(IZ)
      CAI(IZ+1,1)=CAB
      UK(IZ+1)=UKB
      DP(IZ+1,1)=DPB
        DC(IZ+1,1)=DCB
        DB(IZ+1,1)=DBB
        DBT(IZ+1,1)=DBTB
      DO 107 I=1,IZ
        CAI(I,1)=A(5,I)
        CAI(I,2)=A(6,I)
      A(11,I)=ID1(17)
        B(11,I)=ID1(18)
  107 CONTINUE
      IF(TAY.LT..4.AND.K10.GT.1)GOTO 147
   65 YD(1,1)=ID1(14)
      IF(IZ.EQ.1)GOTO 101
   67 YD(IZ,1)=ID1(15)
  101 TE(1)=T1
        TE(IZ+1)=TB
      HB=(DPB-DBB)/2.
        HB1=(DP(1,1)-DB(1,1))/2.
      CALL KHXTA(GKA,HK1,DKH,HKMI,    AKG,IZ)
      IF(K10.GT.1)GOTO 116
  120 IF(ABS(ID1(4)-90.).LT..001)GOTO 116    !判断第一级动叶进口气流角
      CBHA=ID1(27)       !入口导叶相对最大厚度
      XFA=ID1(16)        !入口导叶最大挠度弦相对坐标
        BBA=ID1(19)      !入口导叶迎风度
      BTHA=ID1(21)       !入口导叶稠度
        HBHA=ID1(20)     !入口导叶展弦比
      CALL BXHA(DP1,R1,XFA,BBA,BTHA,HBHA,ALF1*PA,SKB,BHAL,BHAZ)
    6 IF(KDB.EQ.4)GOTO 127         
      GOTO 116
  127 CONTINUE
  116 PI=P1
        TI=T1
        PCY=1.
        A(18,1)=P1
  139 DO 140 I=1,IZ
  142 DP1=DP(I,1)
        DC1=DC(I,1)
        DB1=DB(I,1)
      TI=TE(I)
        C1U=CU(I,1)
        UC1=UK(I)*R0(I,1)
        PI=A(18,I)
      CALL PKIHA(I,UC1,PI,TI,C1U,IZ,PBHA,CMKA)     !转子气动计算
      IF(GB(10).GT.1.)GOTO 147
      PCY=PCY*PC(I,2)
      A(18,I+1)=PI*PC(I,2)
      HZI(I)=A(20,I)/UK(I)**2
      HAD(I)=HZI(I)*XT(I,2)
      A(17,I)=L(I,1)
        B(17,I)=L(I,2)
      HTI(I)=HAD(I)/XT(I,2)/GKA(3,I)
  140 CONTINUE
      DC1=DC(1,1)
        DB1=DB(1,1)
        DP1=DP(1,1)
        UK1=UK(1)
      PCYM=PCY*UBX*PBHA
      XTCYM=(XHA(TXP(XHP(T1)+R*LOG(PCYM)))
     *-XA)/(XHA(TE(IZ+1))-XA)
      UP=ABS(P1 *PCY -PB)/P1 /PCY
      IF(IZ.EQ.1)GOTO147
      IF(IMODD(1).EQ.1)GOTO147
      IF(UP.LT..001) GOTO 147
      K10=2
        XTA=XTCYM
        GOTO 63
  147 IZ7=IZ+1
      DO 2 I=1,IZ7
      DO 2 J=1,2
      AA(I,J)=CAI(I,J)
      AA(I,J+2)=CU(I,J)
      AA(I,5)=AL1(I)
      AA(I,6)=BE1(I)
        AA(I,7)=BE2(I)
      AA(I,8)=AL2(I)
        AA(I,J+8)=AI(I,J)
      AA(I,11)=A(18,I)
        AA(I,12)=UK(I)
        AA(I,13)=HTI(I)
      AA(I,J+13)=L(I,J)
        AA(I,16)=HZI(I)
        AA(I,17)=HAD(I)
      AA(I,18)=GKA(3,I)
        AA(I,19)=A(9,I)
        AA(I,J+19)=PC(I,J)
      AA(I,J+21)=XT(I,J)
        AA(I,J+23)=DL(I,J)
        AA(I,26)=TE(I)
      AA(I,J+49)=B(J,I)
        AA(I,52)=B(15,I)
      AA(I,J+39)=B(J+4,I)
        AA(I,42)=A(12,I)
      AA(I,J+26)=DP(I,J)
        AA(I,J+28)=DB(I,J)
      AA(I,J+30)=DBT(I,J)
        AA(I,J+32)=R0(I,J)
      AA(I,J+34)=F0(I,J)
    2 CONTINUE
      IF(GB(10).GT.1..OR.TAY.LT..4)GOTO 1003
      XS=PCYM/UBX
        X=(XHA(TXP(XHP(T1)+R*ALOG(XS)))-XA)/
     *(XHA(TE(IZ+1))-XA)
  146 IF(K.GT.0.OR.IZ1.GT.0)GOTO 151
      IF(XTCYM.GE.(XTA3-.01 ))GOTO 151
      IF(K11.EQ.0)GOTO 152
  149 IF(   (XTCYM-XTZ).LT..01 )GOTO 151
  152 K11=1
      XTZ=XTCYM
        K10=2
        IDZ=IDZ+1
        IZ=IZ+1
        GOTO 63
  151 CALL HIYX(IZ,KDB)
      LBYX=CAB/SQRT(2.*(XHA(TB)-XHA(TKP(TB))))/SIN(ID1(32)/PA)
      CALL PROH( IZ,A,DP,DB,CAI,CU,UK,CIG,PG,PW,R0)
  145 DO 160 I=1,IZ
      I1=I+1
      A(7,I)=BT(I,1)        !工作轮叶栅稠度
      A(8,I)=BKA(I,1)       !工作轮叶栅长度
      XFK=A(11,I)           !工作轮最大挠度相对坐标
        XFA=B(11,I)         !导向器最大挠度相对坐标
        L1=A(17,I)
        L2=B(17,I)
      KC3=A(10,I)
      DBX=A(14,I)-A(13,I)
        BTK=BT(I,1)
      DAX=B(14,I)-B(13,I)
      CMA=B(12,I)        !导向器叶型面最大相对厚度
        CMK=A(12,I)      !工作轮叶型面最大相对厚度
        BTA=BT(I,2)
        AL4T=B(19,I)*PA
      B2T=B(15,I)
  164 IF(CATAK.EQ.0.)GOTO 162   !攻角和落后角计算
      ZKK=IZKA(I,1)
      ZKA=IZKA(I,2)
      Z11=1.
      CALL ATAKC (DBT(I,1),UK(I),BTK,XFK,CMK,CAI(I,1),
     *CAI(I,2),L1,2.,AL1(I)/PA,BE1(I),BE2(I),B2T,A(5,I),
     *A(6,I),XPC(I,1),AGA1(I,1), AI(I,1),DELT(I,1),TET,
     *ZKK,Z11)                !Howell方法
      UK2=UK(I)*DP(I,2)/DP(I,1)
      CALL ATAKC (DBT(I,2),UK2,BTA,XFA,CMA,CAI(I,2),
     *CAI(I1,1),L2,3.,BE2(I)/PA,AL2(I),AL4(I),AL4T,B(5,I),
     *B(6,I),XPC(I,2),AGA1(I,2),AI(I,2),DELT(I,2),TETA,
     *ZKA,Z11)
      CAI(I,1)=TET
      CAI(I,2)=TETA
      GOTO 163
  162 ZJJ=I
      ZJK=IZ
      DKPI=3.*(AKP-1.)*(ZJJ-.5*ZJK)
      CALL CAB5(KF,KC3,XFK,BE1(I),BE2(I),L1,
     *CMK,BTK,DBX,B2T,DKPI,AI(I,1),DELT(I,1),
     *XPC(I,1),0,AGA1(I,1),CAI(I,1),A(5,I),A(6,I))  !NASA方法
      KC3=0.
      J=2
      CALL CAB5(Kf,KC3,XFA,AL2(I),AL4(I),L2,
     *CMA,BTA,DAX,AL4(I),DKPI,AI(I,2),DELT(I,2),
     *XPC(I,2),2,AGA1(I,2),CAI(I,2),B(5,I),B(6,I))
  163 A(15,I)=0.
        B(15,I)=0.
      CU(I,2)=B(12,I)
        CU(I,1)=A(12,I)
      AA(I,55)=ALX(I)
      DO 160 J=1,2
      AA(I,J+52)=AI(I,J)
  160 CONTINUE
      BHL=BHAL
      DO 166 I=1,IZ
        DO 166 J=1,2
      CAI(I,J)=CAI(I,J)/PA
      L(I,J)=A(J,I)
      B(4,I)=DBT(I+1,1)
        B(7,I)=BT(I,2)
      B(1,I)=DP(I,2)
        B(2,I)=DP(I+1,1)
      A(J,I)=DP(I,J)
        A(3,I)=DBT(I,1)
      DL(I,J)=B(J,I)
      AA(I,J+36)=IZKA(I,J)
      AA(I,39)=A(19,I)
      AA(I,43)=B(12,I)
        AA(I,J+43)=A(J+12,I)
      AA(I,J+45)=B(J+12,I)
        AA(I,J+47)=BKA(I,J)
      AA(I,J+55)=DElT(I,J)
        AA(I,J+57)=AGA1(I,J)
      AA(I,J+59)=YD(I,J)
        AA(I,J+61)=BT(I,J)
      AA(I,J+63)=XPC(I,J)
      I1=I+1
      IF(D3(8).LT..999.AND.A(17,I).GT..85)A(10,I)=1.
      IF(A(17,I).LT..85)GOTO 165
      AI(I,1)=GB(I)
      AI(I,2)=1.
      GOTO 166
  165 AI(I,J)=1.
  166 CONTINUE
      DLKOMP=DB(1,1)
      CALL DLIHA(BKA,AI,DBT,CAI,SKB,BHAL,DLKOMP,DBMAX,KDB,IZ,DS)
      CYMP(20)=UBX
      CYMP(1)=PCYM
        CYMP(2)=XTCYM
        CYMP(3)=DLKOMP
      CYMP(4)=LBYX
        CYMP(5)=HB
        CYMP(6)=HB1
        CYMP(7)=X
      UK1=UK(1)
      PPIB=SQRT(288/T1)
        UK1PP=UK(1)*PPIB
      UKCP=(UK(1)+UK(IZ  ))*.5
        UKCPP=UKCP*PPIB
      GB(1)=P1
        GB(2)=T1
        GB(3)=UK1
        GB(4)=ID1(4)
        GB(6)=R
      CYMP(8)=UK1PP
        CYMP(9)=UKCP
        CYMP(10)=UKCPP
              CYMP(11)=BHAZ
        CYMP(12)=BHL
      GOTO 1002
 1003 IF(GB(10).GT.1.)PRINT 202,GB(9)
 1002 CONTINUE
  201 FORMAT (10X,2HI=,I2,5X,3HF1=,F6.4,5X,3HFB=,F6.4,5X,
     *4HC1A=,F6.2,5X,4HCAB=,F6.2,5X,4HDP1=,F6.4,5X,4HDPB=,F6.4)
  202 FORMAT(10X,24H攢姃帎 剤敂搰帎崕憭?悐 ,F4.3)
      RETURN
        END
      SUBROUTINE ARCH1(SI,DI,AA,IZ,CYMP,OB)
      DIMENSION SI(30),DI(30,2),CYMP(20),AA(30,65),YP(123),X(123),
     *YV(123),q_co(30),fin_speed(4,50,15),fout_speed(4,50,15),
     *twist_const(50,5),fddd(50),det_hu(50),middle_u(50,5),
     *blade_r_in(65,20),blade_r_out(65,20)
      character*100 tempchar
      IZ7=IZ+1
      DO 2 I=1,IZ
      AA(I,17)=AA(I,21)/AA(I,20)
      AA(I,10)=AA(I,1)/AA(I,12)
      AA(I,44)=AA(I,6)+AA(I,53)
      AA(I,45)=AA(I,44)+AA(I,64)
      AA(I,46)=AA(I,8)+AA(I,54)
      AA(I,47)=AA(I,46)+AA(I,65)
      q_co(i)=aa(i,1)/(ob/60*3.1415926*(aa(i,27)+aa(i,29))/2)
    2 continue
c多数参数    

cfdd
      OPEN(7001,FILE='fdd_m1pic.tec') 
      write(7001,*) 'VARIABLES = "stage","fdd"'
      write(7001,*) 'zone t="fdd"'
      do i=1,iz
      write(7001,*) i,aa(i,19)
      enddo
      close(7001)
cc1a,c2a
      OPEN(7012,FILE='c_a.out') 
      do i=1,iz
      write(7012,*) i,aa(i,1),aa(i,2)
      enddo
      close(7012)  
      
      
      
cq_co
      OPEN(7001,FILE='qco_m1pic.tec') 
      write(7001,*) 'VARIABLES = "stage","qco"'
      write(7001,*) 'zone t="qco"'
      do i=1,iz
      write(7001,*) i,q_co(i)
      enddo
      close(7001)   

ctwist
C      OPEN(7003,FILE='s2_in_pre.s12s2')   
	OPEN(7003,FILE='Par_RCu1.1DTOS2')   ! 等环量和等反动度规律的每一列动叶栅进口沿径向12个截面的气流周向速度
      OPEN(7004,FILE='Par_RS12.1DTOS2')      ! 输出动叶栅12个截面的绝对高度位置，用于扭曲设计
	OPEN(7005,FILE='Par_RCU1R.1DTOS2')  ! 每级动叶进口中径处的半径和进口气流周向速度
	OPEN(7009,file='Par_CHORD.1DTOPF') ! 输出中径处弦长用于叶片造型
      do i=1,iz   
      write(7009,*) aa(i,48)
      write(7009,*) aa(i,49)
        rm_in=sqrt((aa(i,27)**2+aa(i,29)**2)/8)
        twist_const(i,1)=aa(i,3)*rm_in
c        twist_const(i,2)=aa(i,4)*(aa(i,28)+aa(i,30))/4
        fddd(i)=aa(i,19)
        middle_u(i,1)=ob/60*3.1415926*2*rm_in
        det_hu(i)=middle_u(i,1)*(aa(i,4)-aa(i,3))
        do j=1,12
cdhl        
            blade_r_in(i,j)=(dble(j-1)/11*(aa(i,27)-aa(i,29))+aa(i,29))
     */2
            blade_r_out(i,j)=(dble(j-1)/11*(aa(i,28)-aa(i,30))+aa(i,30))
     */2

            fin_speed(1,i,j)=twist_const(i,1)/blade_r_in(i,j)
c            fout_speed(1,i,j)=twist_const(i,2)/blade_r_out(i,j)
cdfdd     
            temp_u=ob/60*3.1415926*blade_r_in(i,j)*2
            fin_speed(2,i,j)=temp_u*(1-fddd(i))-det_hu(i)/(2*temp_u)
            temp_u=ob/60*3.1415926*blade_r_out(i,j)*2
            fout_speed(2,i,j)=temp_u*(1-fddd(i))+det_hu(i)/(2*temp_u)
c----------------------------comment after user---------------------     
            
              det_huu=temp_u*(fout_speed(2,i,j)-fin_speed(2,i,j))
              fdddd=1-fin_speed(2,i,j)/temp_u-(fout_speed(2,i,j)-
     !        fin_speed(2,i,j))/2/temp_u

c----------------------------comment after user--------------------- 

        enddo
        
c        rm_out=(aa(i,28)+aa(i,30))/4
        write(7004,'(12f15.8)') (blade_r_in(i,j),j=1,12)
c        write(7004,'(12f15.8)') (blade_r_out(i,j),j=1,12)
        write(7005,'(4f15.8)') aa(i,3),rm_in
      enddo
      

      
      
      do k=1,2
        do i=1,iz
            if(i.lt.10) then
                write(tempchar,'(i1)') i
c               
            else
                write(tempchar,'(i2)') i
            endif
            tempchar=trim(tempchar)//'级'
            
            if(k.eq.1) then
                write(7003,351) '等环量动叶进口第'//tempchar,
     *(fin_speed(k,i,j),j=1,12)
            endif
            if(k.eq.2) then
                write(7003,351) '等反动度动叶进口第'//tempchar,
     *(fin_speed(k,i,j),j=1,12)
            endif
        enddo
      enddo
      
      close(7003) 
      close(7004)
      close(7005)
   
      WRITE(7,301)
      WRITE(7,302)
      WRITE(7,304)
      WRITE(7,1305)
      WRITE(7,310)
      WRITE(7,1305)
      WRITE(7,304)
      WRITE(7,312) (I,(AA(I,J),J=1,13),I=1,IZ7)
      WRITE(7,304)
      WRITE(7,304)
      WRITE(7,305)
      WRITE(7,306)
      WRITE(7,305)
      WRITE(7,304)
      WRITE(7,313) (I,(AA(I,J),J=14,26),I=1,IZ7)
      OPEN(7007,FILE='Par_REFF.1DTOS2')
      OPEN(7008,FILE='Par_TPLS.1DTOS2')   
      
      
      write(7006,*)'PPB'   ! 标识-2
      do iiii=1,IZ7-1
      write(7006,'(3f15.8)') AA(iiii,20),AA(iiii,37),AA(iiii,38)  ! 输出级压比，动叶叶片数，静叶叶片数

        if(iiii.lt.10) then
            write(tempchar,'(i1)') iiii
c               
        else
            write(tempchar,'(i2)') iiii
        endif
        tempchar=trim(tempchar)//'级从叶根到叶顶'
        
        write(7007,351) '动叶等熵效率第'//tempchar,           ! 沿叶高等效率作为S2反设计的初值
     *(AA(iiii,22),j=1,12)
        write(7008,351) '静叶总压恢复系数第'//tempchar,       ! 沿叶高等总压恢复系数作为S2反设计的初值
     *(AA(iiii,17),j=1,12)
      enddo
      close(7007)
      close(7008)
      
      WRITE(7,304)
      WRITE(7,304)
      WRITE(7,303)
      WRITE(7,302)
      WRITE(7,304)
      WRITE(7,305)
      WRITE(7,307)
      WRITE(7,305)
      WRITE(7,304)
      WRITE(7,314) (I,(AA(I,J),J=27,39),I=1,IZ7)
      
      write(7006,*)'HIGHT'     ! 标识-3
      do iiii=1,IZ7           ! 输出各列叶片排出口S2计算站相对位置，沿径向取12个计算栈（可以设置成变量）
        if(iiii.eq.IZ7) then
            write(7006,'(12f15.8)') (dble(jjjj)*         !  dble函数用于转换成双精度
     *((AA(iiii,27)-AA(iiii,29))/22.0),jjjj=0,11)        !  分成11等分，12个计算栈
        else
            write(7006,'(12f15.8)') (dble(jjjj)*         !  动叶片排和静叶片排分开输出
     *((AA(iiii,27)-AA(iiii,29))/22.0),jjjj=0,11)
            write(7006,'(12f15.8)') (dble(jjjj)*
     *((AA(iiii,28)-AA(iiii,30))/22.0),jjjj=0,11)
        endif
      enddo
      
      WRITE(7,304)
      WRITE(7,304)
      WRITE(7,322)
      WRITE(7,305)
      WRITE(7,304)
      WRITE(7,321) (I,(AA(I,J),J=40,52),I=1,IZ7)
      WRITE(7,304)
      WRITE(7,304)
      WRITE(7,305)
      WRITE(7,308)
      WRITE(7,305)
      WRITE(7,304)
      WRITE(7,315) (I,(AA(I,J),J=53,65),I=1,IZ7)
      WRITE(7,304)
      WRITE(7,317) IZ
      WRITE(7,316) (CYMP(I),I=1,7)
      WRITE(7,327)
      WRITE(7,8) (CYMP(I),I= 8,10)
      WRITE(7,1002) (CYMP(I),I=11,12)
      write(7009,*)CYMP(12)
      write(7006,*) 'Num IGV'    ! 标识-4
      write(7006,*) cymp(11)     ! 进口导叶叶片数
      X(1)=0.
      X(2)=CYMP(12)*.75
      X(3)=CYMP(12)
      YV(1)=AA(1,29)/2.
      YV(2)=YV(1)
      YP(1)=AA(1,27)/2.
      YP(2)=YP(1)
      DO 200 J=1,IZ
      X(J*4)=DI(J,1)-SI(J)+X(J*4-1)
      X(J*4+1)=DI(J,1)+X(J*4-1)
      X(J*4+2)=DI(J,2)-SI(J)+X(J*4+1)
      IF(J.EQ.IZ) X(J*4+2)=DI(J,2)+X(J*4+1)
      X(J*4+3)=DI(J,2)+X(J*4+1)
      YV(J*4-1)=AA(J,29)/2.
      YV(J*4)=AA(J,30)/2.
      YV(J*4+1)=YV(J*4)
      YV(J*4+2)=AA(J+1,29)/2.
      YP(J*4-1)=AA(J,27)/2.
      YP(J*4)=AA(J,28)/2.
      YP(J*4+1)=YP(J*4)
      YP(J*4+2)=AA(J+1,27)/2.
  200 CONTINUE
      WRITE(7,201)
  201 FORMAT(5X,'***** 压气机子午流道数据 X,YV,YP *****')
c子午流道
      WRITE(7,202)(X(I),I=1,IZ*4+2)
      WRITE(7,202)(YV(I),I=1,IZ*4+2)
      WRITE(7,202)(YP(I),I=1,IZ*4+2)
      
c子午流道
      OPEN(7001,FILE='zwld_m12s2.tec')
      write(7001,*) 'VARIABLES = "x","yv","yp"'
      write(7001,*) 'zone t="zwld"'
      do i=1,IZ*4+2
      write(7001,*) x(i),yv(i),yp(i)
      enddo
      close(7009)
      close(7001)
      
  351 format(A40,12f15.8)
  202 FORMAT(12F6.4)
  316 FORMAT( 5X,3HPK=,F8.4, 5X,4HKPD=,F6.4, 5X,
     *6HDLINA=,F8.4, 5X,5HLBYX=,F6.4,
     *3X,5HHBYX=,F6.4,3X,4HHBX=,F6.4,
     *3X,6HKPDL4=,F6.4)
  317 FORMAT(20X,19H压气机设计点性能      ,
     *'Z=',I2,1X,10H压气机级数)
  302 FORMAT(40X,38(1H-))
  303 FORMAT(50X,'压气机几何参数列表')
  301 FORMAT(50X,'压气机气动参数列表')
  304 FORMAT(109(1H-))
  305 FORMAT(1HI,3X,13(1HI,7X),1HI)
  306 FORMAT(1HI,1X,1HN,1X,1HI,2X,3HL1T,2X,1HI,3X,2HL3,2X,
     *1HI,3X,2HHZ,2X,1HI,2X,3HSGA,2X,1HI,3X,2HKH,2X,
     *1HI,2X,3HTAY,2X,1HI,3X,2HPK,2X,1HI,3X,2HPC,2X,
     *1HI,1X,4HKPDK,2X,1HI,1X,4HKPDS,2X,1HI,2X,3HDLK,2X,
     *1HI,2X,3HDLA,2X,1HI,3X,2HT1,2X,1HI)
  307 FORMAT(1HI,1X,1HN,1X,1HI,2X,3HDK1,2X,1HI,2X,3HDK2,2X,
     *1HI,2X,3HDB1,2X,1HI,2X,3HDB2,2X,1HI,3X,2HD1,2X,
     *1HI,3X,2HD2,2X,1HI,3X,2HR1,2X,1HI,3X,2HR2,2X,
     *1HI,3X,2HF1,2X,1HI,3X,2HF2,2X,1HI,3X,2HZK,2X,
     *1HI,3X,2HZA,2X,1HI,
     *2X,3HKBT,2X,1HI)
  308 FORMAT(1HI,1X,1HN,1X,1HI,3X,2HIK,2X,1HI,3X,2HIA,2X,
     *1HI,1X,5HALFEK,1X,
     *1HI,1X,5HDELTK,1X,1HI,1X,5HDELTA,1X,1HI,1X,5HAGA1K,1X,
     *1HI,1X,5HAGA1A,1X,1HI,2X,3HYDK,2X,1HI,2X,3HYDA,2X,
     *1HI,2X,3HBTK,2X,1HI,2X,3HBTA,2X,1HI,2X,4HEPSK,1X,
     *1HI,2X,4HEPSA,1X,1HI)
 1305 FORMAT(1HI,3X,8(1HI,7X),1HI,9X,1HI,5X,1HI,9X,2(1HI,7X),1HI)
  310 FORMAT(1HI,1X,1HN,1X,1HI,2X,3HC1A,2X,1HI,2X,3HC2A,2X,
     *1HI,2X,3HC1U,2X,1HI,2X,3HC2U,2X,1HI,2X,3HAL1,2X,
     *1HI,2X,3HBE1,2X,1HI,2X,3HBE2,2X,1HI,2X,3HAL3,2X,
     *1HI,1X,5HENTA1,3X,1HI,1X,3HCAO,1X,1HI,3X,2HP1,4X,1HI,
     *3X,3HUK1,1X,1HI,2X,2HHT,3X,1HI)
  312 FORMAT((1HI,I2,1X,1HI,8(1X,F5.1,2X),F9.0,1X,F5.3,1X,F9.0,1X,
     *F6.1,2X,F6.3,1X,1HI))
  313 FORMAT((1HI,I2,1X,1HI,12(1X,F6.4,1X),F7.2,1HI))
  314 FORMAT ((1HI,I2,1X,1HI,10(1X,F6.4,1X),2(2X,F4.0,2X)
     *,1X,F6.4,1HI))
  315 FORMAT((1HI,I2,1X,1HI,5(F6.2,2X),8(F7.4,1X),1HI))
  321 FORMAT((1HI,I2,1X,1HI,4(1X,F6.4,1X),4(1X,F6.2,1X),
     *2(1X,F6.4,1X),2(F8.0,1X),F6.2,1X,1HI))
  322 FORMAT(1HI,1X,1HN,1X,1HI,1X,4HZITK,2X,1HI,
     *1X,4HZITA,2X,1HI,2X,2HCK,
     *3X,1HI,2X,2HCA,3X,1HI,2X,3HB1E,2X,1HI,2X,3HB2E,2X,
     *1HI,1X,4HAL3E,2X,1HI,1X,4HAL4E,2X,1HI,2X,2HBK,3X,
     *1HI,2X,2HBA,3X,1HI,2X,3HP1C,2X,1HI,2X,3HP2C,2X,1HI,
     *2X,3HB2T,2X,1HI)
  327 FORMAT(5X,6HUK1 PR,9X,4HUKCP,5X,7HUKCP PR)
 1002 FORMAT(5X,5HZBHA=,F3.0,5X,5HLBHA=,F6.4)
 1314 FORMAT(1X,1HI,
     * 5X,4HDLPK,10X,4HDLHA,5X,5H噣噹?/,
     *(I2,3(F10.5,4X)))
 1315 FORMAT(5X,18H姁垝厫垑 彁帡崕憭?/,1X,1HI,1 X,3HKP1,
     *5 X,3HKP2,5 X,3HKP3,/,(I2,3(I2,5 X)))
   50 FORMAT(5X,10H寑憫€ 弲?,F8.1/5X,10H寑憫€ 妿?,F8.1,
     *5X,10H寑憫€ 憮?,F8.1)
  326 FORMAT(5X,10H寑憫� 倣�=,F8.1)
   18 FORMAT(5X,6HDLBHA=,F7.4,5X,6HDLPER=,F7.4)
    8 FORMAT(5F13.4)
 1316 FORMAT (5X,10HSIGMA PER=,F6.4)
      RETURN
      END
      SUBROUTINE PERKSI(AA,BXD1,IZ,PKSI,DLP,RV,DLC,GAM,ALP,ALF1)
      DIMENSION BXD1(10),AA(30,65)
      REAL N
      IZZ=IZ+1
      N=BXD1(2)
      DWT0=AA(IZZ,29)
      DWT1=BXD1(7)
      H0=(AA(IZZ,27)-DWT0)/2.
      ALF0=AA(IZZ,5)/57.296
      B0=(DWT0+H0)/N*3.1416*SIN(ALF0)
      CTR=.0300
      F0=3.1416*(AA(IZZ,27)**2-DWT0**2)*SIN(ALF0)/N/4.
      ALP=BXD1(3)/57.296
      IF(ALP.EQ.0) ALF1=ALF0
      GAM=ATAN((AA(IZ,30)-DWT0)/2./AA(IZ,49))
      IF(ALP.EQ.0.) GOTO 6
      IF(GAM.GT.ALP) ALP=GAM
      IF(DWT1.LT.DWT0) GOTO 1
      PRINT*,'DWT1.GE.DWT0, DWT0=',DWT0
      STOP
    1 CONTINUE
      DR=(DWT0-DWT1)/2.
      RV=DR/(COS(GAM)-2.*COS(ALP)+1.)
      DL=RV*(2.*ALP-GAM)
      GOTO 7
    6 CONTINUE
      IF(ABS(GAM).LT..0173) GAM=.0173
      DL=AA(IZ,49)*.25
      RV=ABS(DL/SIN(GAM))
      IF(GAM.LT.0) DWT1=DWT0+2.*RV*(1.-COS(GAM))
      IF(GAM.GT.0) DWT1=DWT0-2.*RV*(1.-COS(GAM))
    7 CONTINUE
      RS=RV+H0/2.
      DL=DL*RS/RV
      R0=SQRT(F0/3.1416)
      FIKR=0
      DO 3 I=1,100
      FIKR=FIKR+.00873
      R1=R0+DL*TAN(FIKR/2.)
      F1=3.1416*R1**2
      FIK=.43*(CTR*(F1+F0)/(F1-F0))**.4444
      IF(FIK.LE.FIKR) GOTO 4
    3 CONTINUE
    4 CONTINUE
      PSR1=3.1416*(R0+R1)
      FISR=2.*(F1-F0)/DL/PSR1
      H1=(SQRT((AA(IZZ,27)**2-DWT0**2)*SIN(ALF0)/SIN(ALF1)+
     *DWT1**2)-DWT1)*.5
      DP1=SQRT(F0*4./3.1416/SIN(ALF1)*N+DWT1**2)
      H1=(DP1-DWT1)/2.
      P0=2.*(H0+B0)
      DO 30 I=1,1000
      DP1=DP1+.01*H0
      H1=(DP1-DWT1)/2.
      F1=3.1416/4.*(DP1**2-DWT1**2)*SIN(ALF1)/N
      B1=(DWT1+H1)/N*3.1416*SIN(ALF1)
      P1=2.*(H1+B1)
      PSR2=(P0+P1)/2.
      FIS=2.*(F1-F0)/DL/PSR2
      IF(FIS.GE.FISR) GOTO 40
   30 CONTINUE
   40 CONTINUE
      IF(ALP.EQ.0.) DP1=SQRT(F0*4./3.1416/SIN(ALF1)*N+DWT1**2)
      GOTO 24
   25 PRINT*,'搨厠垪湌?DP1'
   24 PRINT*,'DP1=',DP1
      PRINT*,'DP1='
      READ(5,*) DP1
      BXD1(6)=DP1
      BXD1(7)=DWT1
      BXD1(5)=BXD1(7)/BXD1(6)
      H1=(DP1-DWT1)/2.
      B1=(DWT1+H1)/N*3.1416*SIN(ALF1)
      F1=3.1416/4.*(DP1**2-DWT1**2)*SIN(ALF1)/N
      HSR=(H0+H1)/2.
      BSR=(B0+B1)/2.
      DG=2.*HSR*BSR/(HSR+BSR)
      IF(ALP.EQ.0) GOTO 8
      A11=.9*SIN(ALP-GAM)
      A12=.9*SIN(ALP)
      A1=(A12+A11)*(1.+.5*(1.-GAM/ALP))
      GOTO 9
    8 CONTINUE
      A1=.9*SIN(GAM)
    9 CONTINUE
      B1=.21/SQRT(RS/DG)
      BH=BSR/HSR
      IF(BH.LT..5.OR.BH.GT.8.) PRINT 2,BH
    2 FORMAT(' WARNING!  BH=',F5.2)
      IF(BH.LE.3.) C1=(3.-BH)**3.585*.05+.4
      IF(BH.GT.3.) C1=.1025*SIN(.5712*BH-2.9988)+.4983
      ZM=A1*B1*C1
      IF(GAM.EQ..0173) ZM=0.
      ZT=.0175*CTR*RS/DG*(2.*ALP-GAM)*57.296
      IF(ALP.EQ.0.) ZT=-ZT
      TG=H0*(F1-F0)/F0/2./DL
      ZR=3.2*(ABS(TG))**1.25*(1.-F0/F1)**2
      ZS=ZT+ZR+ZM
      QL0=1.5774*(1.-.1667*AA(IZ,15)**2)**2.5*AA(IZ,15)
      QL1=QL0*F0/F1
      IF(QL1.GE.1.) GOTO 25
      AL1=ALQ(QL1,1.4)
      ALS=(AA(IZ,15)+AL1)/2.
      PKSI=1.-.5833*ALS**2*ZS
      DLP=RV*(2.*SIN(ALP)-SIN(GAM))*SIN(ALF0)
      IF(ALP.EQ.0) DLP=-DLP
      IF(DLP.LT.(AA(IZ,49)*.25)) DLP=AA(IZ,49)*.25
      DLC=RV*(.5-SIN(GAM))*SIN(ALF0)
      WRITE(7,11)
   11 FORMAT(5X,'********** 悁憲厭 弲悈晭剭垔€ **********')
      GA=GAM*57.3
      WRITE(7,12)GA,RV,DLP,PKSI
      WRITE(7,13)
   13 FORMAT(5X,'妿帎剤崁挍 弲悈晭剭垔€  X,YV,YP')
   12 FORMAT(5X,4HGAM=,F5.1,5X,3HRV=,F6.3,5X,4HDLP=,F5.3,
     *5X,5HPKSI=,F5.4)
      RETURN
      END
      FUNCTION ALQ(QL,AK)
      X=.6*QL
      Y=1./(AK-1.)
    3 CONTINUE
      RLI=QL/(((AK+1.-(AK-1.)*X**2)/2.)**Y)
      Z=ABS(RLI-X)
      IF(Z-1E-6)2,2,1
    1 X=RLI
      GOTO 3
    2 ALQ=RLI
      RETURN
      END
      SUBROUTINE PER(AA,BXD1,BXD2,IZ,DLP,G3,RV,PIK,GAM,ETAS,SI,DI,SYMP,
     *PKSI,ALP,ETCBK,A1C,ALF1,IOCK,OB,PSI,GAMLD)
      DIMENSION SI(30),DI(30,2),AA(30,60),BXD1(10),BXD2(50),SYMP(20),
     *S2(20),S1(20),X2(50),Y2(50),Y3(50),X1(11),Y1(11),X3(10)
      DATA X1/0.,.14,.28,.41,.54,.65,.75,.84,.91,.96,.99/,
     *Y1/0.,.00985,.04,.0879,.1583,.24,.3386,.4574,.5854,.72,.8589/
      PSI=PSI/57.296
      IG=47
      IB=1
  101 CONTINUE
      DLCD=BXD2(6)+BXD2(6)*TAN(PSI)
      DOCK=SYMP(3)+DLP+DLCD
      IF(IOCK.EQ.2) GOTO 100
      DO 22 L=1,IZ
      S2(L+1)=SI(L)
      S1(L)=SI(L)
   22 CONTINUE
      X=SYMP(12)
      S2(1)=X*.25
      X2(1)=0.
      Y2(1)=AA(1,27)*.5
      X2(2)=X*.75
      Y2(2)=Y2(1)
      X2(3)=X2(2)
      Y2(3)=AA(1,29)*.5
      X2(4)=0.
      Y2(4)=Y2(3)
      X=SYMP(12)-S2(1)
      Y=AA(1,29)*.5
      DO 2 N=1,IZ
      X2(1)=X+S2(N)
      Y2(1)=Y
      X2(2)=X2(1)
      Y2(2)=AA(N,27)*.5
      X2(3)=X2(1)+DI(N,1)-S1(N)
      Y2(3)=AA(N,28)*.5
      X2(4)=X2(3)
      Y2(4)=0.
      X2(5)=X2(1)
      Y2(5)=0.
      X2(6)=X2(1)
      Y2(6)=Y2(1)
      X2(7)=X2(3)
      Y2(7)=AA(N,30)*.5
      X2(8)=X2(7)+S1(N)+DI(N,2)-S2(N+1)
      IF(N.EQ.IZ) X2(8)=X2(7)+S1(N)+DI(N,2)
      Y2(8)=AA(N+1,29)*.5
      X2(9)=X2(8)
      Y2(9)=AA(N+1,27)*.5
      X2(10)=X2(9)-DI(N,2)+S2(N+1)
      IF(N.EQ.IZ) X2(10)=X2(9)-DI(N,2)
      Y2(10)=AA(N,28)*.5
      X2(11)=X2(10)
      Y2(11)=AA(N,30)*.5
      X=X2(8)
      Y=Y2(8)
    2 CONTINUE
      DLS=X
      IF(IOCK.EQ.1) GOTO 27
      X2(1)=X
      Y2(1)=Y
      Y3(1)=AA(IZ+1,27)*.5
      X2(12)=X2(1)+DLP
      Y2(12)=BXD1(7)*.5
      Y3(12)=BXD1(6)*.5
      ALF0=AA(IZ+1,5)/57.296
      H0=Y3(1)-Y2(1)
      H1=Y3(12)-Y2(12)
      DH=(H1-H0)/11.
      DHS=0.
      IF(ALP.GT.0.) GOTO 17
      TET=ABS(GAM)
      DTET=TET/11.
      DO 18 I=2,11
      TET=TET-DTET
      IF(GAM.LE.0.) Y2(I)=Y2(12)-RV*(1.-COS(TET))
      IF(GAM.GT.0.) Y2(I)=Y2(12)+RV*(1.-COS(TET))
      X2(I)=X2(12)-RV*SIN(TET)*SIN(ALF0)
      DHS=DHS+DH
      Y3(I)=Y2(I)+H0+DHS
   18 CONTINUE
      GOTO 9
   17 CONTINUE
      IF(GAM.GE.0.) GOTO 6
      TET1=0.
      DO 7 I=2,7
      Y2(I)=Y2(12)+RV*2.*(1.-COS(ALP))-RV*(1.-COS(TET1))
      X2(I)=X2(12)-(RV*2.*SIN(ALP)-RV*SIN(TET1))*SIN(ALF0)
      TET1=TET1+.2*ALP
      DHS=DHS+DH
      Y3(I)=Y2(I)+DHS+H0
    7 CONTINUE
      TET2=.8*ALP
      DO 8 I=8,11
      Y2(I)=Y2(12)+RV*(1.-COS(TET2))
      X2(I)=X2(12)-RV*SIN(TET2)*SIN(ALF0)
      TET2=TET2-.2*ALP
      DHS=DHS+DH
      Y3(I)=Y2(I)+DHS+H0
    8 CONTINUE
      GOTO 9
    6 CONTINUE
      N1=5*(1.-GAM/ALP)
      IF(N1.EQ.0) N1=1
      N2=11-N1
      DTET1=(ALP-GAM)/N1
      TET1=GAM+DTET1
      N3=N1+1
      IF(N1.GE.2) GOTO 19
      N2=11
      N3=2
      DTET2=ALP/N2
      TET2=ALP-DTET2
      GOTO 12
   19 CONTINUE
      DO 10 I=2,N1
      Y2(I)=Y2(12)+RV*2.*(1.-COS(ALP))-RV*(1.-COS(TET1))
      X2(I)=X2(12)-(RV*2.*SIN(ALP)-RV*SIN(TET1))*SIN(ALF0)
      TET1=TET1+DTET1
      DHS=DHS+DH
      Y3(I)=Y2(I)+DHS+H0
   10 CONTINUE
      DTET2=ALP/N2
      TET2=ALP
   12 CONTINUE
      DO 11 I=N3,11
      Y2(I)=Y2(12)+RV*(1.-COS(TET2))
      X2(I)=X2(12)-RV*SIN(TET2)*SIN(ALF0)
      TET2=TET2-DTET2
      DHS=DHS+DH
      Y3(I)=Y2(I)+DHS+H0
   11 CONTINUE
    9 CONTINUE
      WRITE(7,37)(X2(I),I=1,12)
      WRITE(7,37)(Y2(I),I=1,12)
      WRITE(7,37)(Y3(I),I=1,12)
   37 FORMAT(12F6.4)
  100 CONTINUE
      IF(IOCK.EQ.2) X2(12)=0.
      SVC=(BXD1(6)-BXD1(7))/4.
      X2(1)=X2(12)
      Y2(1)=BXD1(6)/2.
      X2(2)=X2(1)+SVC*.75
      Y2(2)=Y2(1)
      X2(3)=X2(2)
      Y2(3)=BXD1(7)/2.
      X2(4)=X2(1)
      Y2(4)=Y2(3)
      X2(5)=X2(4)
      Y2(5)=Y2(1)
      IF(A1C.EQ.ALF1) GOTO 20
   20 CONTINUE
      X2(1)=X2(2)+SVC*.25
      IF(A1C.EQ.ALF1) X2(1)=X2(4)
      Y2(1)=Y2(3)
      X2(2)=X2(1)
      Y2(2)=BXD1(6)/2.
      HP=(BXD2(16)-BXD1(6))/2.+BXD2(17)/2.*SIN(PSI)
      BP=BXD2(6)-BXD2(17)*COS(PSI)-HP*TAN(PSI)
      BPP=BP
      IF(BPP.LT.HP*.7) BP=HP*.7
      IF(BPP.LT.HP*.7) BXD2(6)=HP*TAN(PSI)+BP+BXD2(17)*COS(PSI)
      DO 13 I=3,12
      Y2(I)=HP*Y1(I-1)+Y2(2)
      X2(I)=BP*X1(I-1)+X2(2)+(Y2(I)-Y2(2))*TAN(PSI)
   13 CONTINUE
      X2(13)=X2(2)+BP+HP*TAN(PSI)
      Y2(13)=BXD2(16)/2.+BXD2(17)/2.*SIN(PSI)
      X2(14)=X2(13)+BXD2(17)*COS(PSI)
      Y2(14)=Y2(13)-BXD2(17)*SIN(PSI)
      DUG=.157
      IF(PSI.LT.DUG) DUG=PSI
      HV=(BXD2(16)-BXD1(7))/2.-BXD2(17)/2.*SIN(PSI)
      BV=BXD2(6)-HV*TAN(PSI-DUG)
      DO 14 I=15,24
      Y2(I)=Y2(1)+Y1(26-I)*HV
      X2(I)=X2(1)+X1(26-I)*BV+(Y2(I)-Y2(1))*TAN(PSI-DUG)
   14 CONTINUE
      X2(25)=X2(1)
      Y2(25)=Y2(1)
      X2(26)=X2(25)
      Y2(26)=0.
      X2(27)=X2(26)+BXD2(6)
      Y2(27)=0.
      X2(28)=X2(27)
      Y2(28)=Y2(14)
      WRITE(7,38)
   38 FORMAT(5X,'*** 妿帎剤崁挍 枀崚悗亝啀巸?妿媴憖  XP,YP;XV,YV ***')
      WRITE(7,37)(X2(I),I=2,13)
      WRITE(7,37)(Y2(I),I=2,13)
      WRITE(7,37)(X2(I),I=14,25)
      WRITE(7,37)(Y2(I),I=14,25)
      X2(1)=X2(13)
      Y2(1)=Y2(13)
      X2(2)=X2(1)+(BXD2(26)-BXD2(16))/2.*TAN(PSI)
      Y2(2)=Y2(1)+(BXD2(26)-BXD2(16))/2.
      X2(3)=X2(1)+(BXD2(34)-BXD2(16))/2.*TAN(PSI)
      Y2(3)=Y2(1)+(BXD2(34)-BXD2(16))/2.
      Y2(7)=BXD2(40)/2.
      X2(7)=X2(3)+(Y2(7)-Y2(3))*TAN(PSI)+(Y2(7)-Y2(3))/COS(PSI)
      R34=(Y2(7)-Y2(3))/COS(PSI)*TAN(1.5708-(1.5708-PSI)/2.)
      HOR=2.*R34*SIN((1.5708-PSI)/8.)
      DO 15 I=4,6
      X2(I)=X2(I-1)+HOR*SIN(PSI+(1.5708-PSI)/8.*((I-4)*2.+1.))
      Y2(I)=Y2(I-1)+HOR*COS(PSI+(1.5708-PSI)/8.*((I-4)*2.+1.))
   15 CONTINUE
      X3(1)=X2(1)+BXD2(17)*COS(PSI)
      Y3(1)=Y2(1)-BXD2(17)*SIN(PSI)
      X3(2)=X2(2)+BXD2(30)*COS(PSI)
      Y3(2)=Y2(2)-BXD2(30)*SIN(PSI)
      X3(3)=X2(3)+BXD2(35)*COS(PSI)
      Y3(3)=Y2(3)-BXD2(35)*SIN(PSI)
      X3(7)=X2(7)
      Y3(7)=BXD2(41)/2.
      PSV=PSI+GAMLD
      R35=(Y3(7)-Y3(3))/COS(PSV)*TAN(1.5708-(1.5708-PSV)/2.)
      HOR=2.*R35*SIN((1.5708-PSV)/8.)
      DO 16 I=4,6
      X3(I)=X3(I-1)+HOR*SIN(PSV+(1.5708-PSV)/8.*((I-4)*2.+1.))
      Y3(I)=Y3(I-1)+HOR*COS(PSV+(1.5708-PSV)/8.*((I-4)*2.+1.))
   16 CONTINUE
      X2(1)=X2(2)
      Y2(1)=Y2(2)
      X2(2)=X3(2)
      Y2(2)=Y3(2)
      X2(1)=X2(3)
      Y2(1)=Y2(3)
      X2(2)=X3(3)
      Y2(2)=Y3(3)
      DLS=X3(7)
   27 CONTINUE
      BXD1(5)=BXD1(7)/BXD1(6)
      IF(IB.EQ.2) GOTO 102
      IB=1
      IG=32
      IF(IB.EQ.2) GOTO 101
  102 CONTINUE
      RETURN
      END
      SUBROUTINE GH(NV,GPR,PI,ET,TETA,IZ,AN,PK,GPRR,IOCK)
      DIMENSION GPR(10,40),PI(10,40),ET(10,40),X(40),Y(40),TETA(10,30)
     *,AN(10,3)
      GMAX=GPRR
      GMIN=10000.
      DO 10 I=1,10
      DO 10 J=1,40
      IF(GPR(I,J).GT.GMAX) GMAX=GPR(I,J)
      IF(GPR(I,J).GT.0..AND.GPR(I,J).LT.GMIN) GMIN=GPR(I,J)
   10 CONTINUE
      IZ1=IZ+1
      IG=47
      IB=1
  101 CONTINUE
      IF(PK.GT.1.) YMAX=4.
      IF(PK.GT.2.8) YMAX=8.
      IF(PK.GT.6.) YMAX=16.
      IF(PK.GT.12.) YMAX=32.
      UY=YMAX/8.
      UUX=GMAX/8.
      IF(UUX.GT..1) UX=.1
      IF(UUX.GT..2) UX=.2
      IF(UUX.GT..4) UX=.4
      IF(UUX.GT.1.) UX=1.
      IF(UUX.GT.2.) UX=2.
      IMAX=GMAX/UX
      TMAX=IMAX+1
      XMAX=UX*TMAX
      XMIN=XMAX-UX*8.
      DO 1 N=1,NV
      J1=0
      NM=N
      DO 2 J=1,40
      IF(GPR(N,J).EQ.0..OR.PI(N,J).LT.1.) GOTO 2
      J1=J1+1
      X(J1)=GPR(N,J)
      Y(J1)=PI(N,J)
    2 CONTINUE
    1 CONTINUE
      DO 3 N=1,NV
      J2=0
      NM=N
      DO 4 J=1,40
      IF(GPR(N,J).EQ.0..OR.ET(N,J).LT..7.OR.PI(N,J).LT.1.) GOTO 4
      J2=J2+1
      X(J2)=ET(N,J)
      Y(J2)=PI(N,J)
    4 CONTINUE
    3 CONTINUE
      XM=1.0
      YM=16.2
      YMM=16.4
      DO 5 N=1,NV
      NM=N
      YM=YM-1.
      YMM=YMM-1.
      ANM=AN(N,1)
      XM=XM+1.2
      DO 6 I=1,IZ1
      TET=TETA(N,I)
      XM=XM+2.
    6 CONTINUE
      XM=1.0
    5 CONTINUE
      YM=YM-1.
      ZZ=IZ
      IF(IOCK.EQ.2) ZZ=0.
      YM=YM-.8
      YM=YM-.2
      XM=2.2
      IF(IB.EQ.2) GOTO 102
      IB=1
      IG=32
      IF(IB.EQ.2) GOTO 101
  102 CONTINUE
      RETURN
      END
      SUBROUTINE PE4ATI(IPE4,KPAH,FKPAH,KRAN,II,AN,M,P,AKPX,IG)
      DIMENSION KPAH(40),KRAN(40),FKPAH(10),AN(10,3)
      INTEGER M
      FI=II
      IF(IPE4.NE.1) RETURN
      IF(AN(M,1).LT.FKPAH(1).OR.AN(M,1).GT.FKPAH(2)) GOTO 1       ! FKPAH 数值均为0
      IF(AN(M,2).LT.FKPAH(3).OR.AN(M,2).GT.FKPAH(4)) GO TO 1
      IF(IG.EQ.1) GOTO 4
      IF(P.LT.FKPAH(5).OR.P.GT.FKPAH(6)) GO TO 1
      IF(AKPX.LT.FKPAH(7).OR.AKPX.GT.FKPAH(8)) GO TO 1
 4    IF(FI.LT.FKPAH(9).OR.FI.GT.FKPAH(10)) GO TO 1
      DO 2 L=1,40
    2 KRAN(L)=KPAH(L)
      RETURN
  1   DO 3 L=1,40
    3 KRAN(L)=0
      RETURN
      END
      SUBROUTINE STEPON(M,IQP,IQ,IDQC,INDEXQ,INDEX2,INDEX4,INDEX6,INDEX7
     *,IPE4,GB,AN,DELTQ,P,QL1,AKPX,QH,QGP,PGP,QQS,QQF,IM,
     *KM,INDEX8,IND8R,INDEXB,KSTEP,BRANCH,POINTS)
      DIMENSION GB(10),AN(10,3), BRANCH(10)
      INTEGER M
 700  FORMAT(16I5)
 770  FORMAT(10F12.4)
      KSTEP=0
      IF(IQP.EQ.0) GO TO 1        ! IQP = 0
      IF(IQP.EQ.1) GO TO 41
      IF(INDEX2.EQ.0) GO TO 42
      IF(IQP.EQ.3) GO TO 33
      IF(IQ.EQ.0) IQ=1
      AN(M,2)=AN(M,2)-ABS(DELTQ)
      IM=0
      KM=0
      INDEX2=0
      KSTEP=5
      RETURN
   42 IF(IQP.EQ.2.AND.IQ.EQ.1) IQP=3
      IF(IQP.EQ.3) GO TO 45
      KSTEP=4
      RETURN
   45 CONTINUE
      IF(INDEX7.EQ.1.OR.INDEX6.EQ.1) GO TO 33
      INDEX4=1
      IND8R=1
      AN(M,2)=AN(M,2)+ABS(DELTQ)
      GO TO 12
   41 CONTINUE
      IF(INDEX2.EQ.0) GO TO 2
      IF(IPE4.EQ.0) GO TO 3
      WRITE(6,100)
      WRITE(6,101)
  100 FORMAT('    确定阻塞工况点')
  101 FORMAT('    不存在的计算')
    3 POINTS=1
    2 KSTEP=3
      RETURN
    1 IF(IND8R.EQ.1) GO TO 33
      IF(IDQC.NE.0.AND.IDQC.LT.IM*INDEX2) GO TO 4
   50 IF(INDEX4.NE.0.AND.INDEX2.NE.0) GO TO 33
      IF(INDEX7.EQ.0.AND.INDEX6.EQ.0) GO TO 5
   33 QH=QH+(0.5-INDEX2)*ABS(DELTQ)*(1-INDEX6)
      KSTEP=2
      RETURN
    5 IF(INDEX4.NE.0) GO TO 7
      IF(INDEX2.EQ.0) GO TO 8
      IF(AN(M,2).GT.QQS) GO TO 27
      IF(IPE4.EQ.0) GO TO 2
      WRITE(7,104)
      WRITE(7,105)
  104 FORMAT('    给定融冰的流量')
  105 FORMAT('    压气机喘振工作')
      GO TO 2
   27 CONTINUE
      KSTEP=5
      INDEX2=0
      GO TO 9
    8 INDEX4=1
      IF(QQS.LT.1E-6) GO TO 9
      DELTQ=ABS(DELTQ)
      GO TO 11
    9 QH=AN(M,2)
      QGP=QH
      PGP=P
   10 DELTQ=-ABS(DELTQ)
   11 AN(M,2)=AN(M,2)+DELTQ
      IF(AN(M,2).GT.0.) GO TO 26
      INDEXQ=2
      AN(M,2)=AN(M,2)+2.*ABS(DELTQ)
   26 IF(AN(M,2).LT.1) GO TO 25
      KSTEP=4
      RETURN
   25 CONTINUE
      IF(QQF.LT.1E-6) GO TO 20
      IF(AN(M,2).GT.QQF) AN(M,2)=QQF
   20 CONTINUE
   12 IF(KSTEP.NE.5) KSTEP=1
      RETURN
    7 IF(INDEXQ.NE.0) GO TO 21
      IF(QQS.GT.1E-6) GO TO 11
      IF(P.LT.PGP) GO TO 13
   14 PGP=P
      QGP=AN(M,2)
      INDEXQ=1
      GO TO 10
   13 INDEXQ=2
   15 DELTQ=ABS(DELTQ)
      AN(M,2)=QH+DELTQ
      QH=AN(M,2)
      GO TO 12
   21 IF(INDEXQ.EQ.1) GO TO 16
      IF(P.LT.PGP) GO TO 15
      PGP=P
      QGP=AN(M,2)
      GOTO 15
   16 IF(P.GE.PGP.AND.INDEXB.NE.4) GO TO 14
      GO TO 15
    4 IF(INDEX6.EQ.1) GO TO 33
      IF(ABS(DELTQ).GT.GB(9)/10.) GO TO 17
      BRANCH(M)=AN(M,1)
      IF(IPE4.EQ.0) GO TO 2
      WRITE(7,102)
      WRITE(7,103)
  102 FORMAT('    给定的条件不能工作')
  103 FORMAT('    失速与喘振同时存在')
      GO TO 2
   17 INDEX4=1
      INDEXQ=2
      IM=0
      KM=0
      INDEXB=0
      INDEX2=0
      KSTEP=5
      QH=AN(M,2)+ABS(DELTQ)
      IF(INDEX7.NE.0) DELTQ=0.5*ABS(DELTQ)
      DELTQ=ABS(DELTQ)
      GO TO 11
      END
      SUBROUTINE STEPAT(INDEX2,INDEXB,INDEXD,INDEXQ,INDEX6,INDEX8,
     *IND8R,QQS,QQF,E,ES,QL1,P,G,Q1,P1,QGP,PGP,GGP,Q2,P2,DELTQ,
     *QMEM,DQMEM,ZAP,IQP,KSTEP)
      KSTEP=0
      IF(IQP*INDEXB.EQ.12) GOTO 1
      IO=INDEXB/4*((P+1)/(PGP+1))
      INDEXB=INDEXB*(1-IO)
      IF(IQP.GT.0) GOTO 1
      IF(INDEX6.EQ.1) GO TO 1
      IF(QQS.GT.1E-8.AND.QQF.GT.1E-8) GO TO 1
      IF(INDEX2.EQ.1) GO TO 1
      IF(INDEXB.EQ.4) GO TO 1
      IF(INDEXB.GE.2) GO TO 2
      ES=E*0.2
      IF(INDEX8.GE.1.OR.IND8R.EQ.1) ES=E/10.
      IF(QL1.LT.QGP) GO TO 3
      IF(P.LE.PGP) GO TO 4
      Q1=QGP
      P1=PGP
      QGP=QL1
      PGP=P
      GGP=G
      ZAP=P/G
      GO TO 1
    4 QMEM=QL1
      DQMEM=DELTQ
      Q2=QL1
      P2=P
      INDEXB=2
      INDEXD=2
      GO TO 2
    3 IF(P.LE.PGP) GO TO 5
      Q2=QGP
      P2=PGP
      QGP=QL1
      PGP=P
      GGP=G
      ZAP=P/G
      GO TO 1
    5 Q1=QL1
      P1=P
      IF(INDEXQ.EQ.0) GO TO 1
      INDEXB=2
      DQMEM=DELTQ
      QMEM=QL1
      GO TO 2
    1 KSTEP=1
      RETURN
    2 IF(INDEXB.NE.3 ) GO TO 6
      INDEXB=4
      QL1=QMEM
      DELTQ=DQMEM
      GO TO 1
    6 IF(Q2-Q1.LT.ES) INDEXB=3
      IF(INDEXD.NE.0) GO TO 7
      IF(P.LT.PGP) GO TO 8
      Q2=QGP
      P2=PGP
      QGP=QL1
      PGP=P
      GGP=G
      ZAP=P/G
      GO TO 9
    8 Q1=QL1
      P1=P
    9 IF(Q2-QGP.LT.ES/2.) GO TO 10
   14 QL1=(Q2+QGP)/2.
      INDEXD=2
      GO TO 11
    7 IF(P.LT.PGP) GO TO 12
      Q1=QGP
      P1=PGP
      QGP=QL1
      PGP=P
      GGP=G
      ZAP=P/G
      GO TO 13
   12 Q2=QL1
      P2=P
   13 IF(QGP-Q1.LT.ES/2.) GO TO 14
   10 QL1=(Q1+QGP)/2.
      INDEXD=0
   11 KSTEP=2
      RETURN
      END
      SUBROUTINE DECOD(COD,K,A,N,B,M)
        DOUBLE PRECISION COD,AI,COD1,A1
        INTEGER*4 I
        COD1=COD/10**K
        I=COD1
      AI=DFLOAT(I)
      A1=AI/10**N
      A=SNGL(A1)
      AI=AI*10**K
      B=(COD-AI)/10**M
      RETURN
      END
      SUBROUTINE STEPIN(INDEXZ,INDEX2,INDEX6,INDEX9,INDEXB,LSR,IM,KM,IQP
     *,PR,IZAP,CU,ZAP,QL1,G,P,DELTQ,SIGMA,DSIGMA,QSET,DQSET,SISET,DSISET
     *,SDELTQ,WCP,EP,GUIDE,SGUIDE,KSTEP)
      DIMENSION WCP(10,4)
      KSTEP=0
      IF(INDEXZ.NE.0) GO TO 1
	KU=INDEX6+1
      IF(IQP.GE.2.AND.LSR.EQ.0) GOTO (20,10),KU
      IF(IZAP.EQ.0) GO TO 20
      IF(INDEXB.NE.4) GO TO 20
	PZ=G*ZAP/(1.+CU/100.)
      IF(LSR.EQ.1) PZ=GWCP(G,WCP)
   10 IF(IQP.GE.2.AND.LSR.EQ.0) PZ=PR
      IF(P.GT.PZ) GO TO 20
      INDEXZ=1
      IF(INDEX6.EQ.0) GO TO 2
      IM=1000
      KM=2
    2 QSET=QL1
      DQSET=DELTQ
      SISET=SIGMA
      DSISET=DSIGMA
      SGUIDE=GUIDE
      IF(SDELTQ.EQ.0) SDELTQ=DELTQ
      GO TO 30
    1 IF(INDEXZ.GE.2) GO TO 20
    3 PZ=G*ZAP/(1.+CU/100.)
      IF(LSR.EQ.1) PZ=GWCP(G,WCP)
      IF(IQP.GE.2.AND.LSR.EQ.0) PZ=PR
   30 CONTINUE
      IF(INDEX6.EQ.1) GO TO 4
      IF(ABS(P-PZ).LT.EP) GO TO 5
      IF(ABS(SDELTQ).LT.1E-6) GOTO 5
      SDELTQ=ABS(SDELTQ)/2.
      IF(P.LT.PZ) SDELTQ=-2*SDELTQ
      QL1=QL1+SDELTQ
      GO TO 15
    5 INDEXZ=2
      QL1=QSET
      DELTQ=DQSET
      GO TO 20
    4 IF(ABS(P-PZ).LT.EP) GO TO 6
      DSIGMA=ABS(DSIGMA)/2.
      IF(DSIGMA.LE.1.E-6) GOTO 20
      IF(P.LT.PZ) DSIGMA=-DSIGMA
      INDEX9=1
      GO TO 7
    6 INDEXZ=2
      SIGMA=SISET
      DSIGMA=DSISET
      GUIDE=SGUIDE
      IF(IQP.EQ.3) GO TO 21
      IM=0
      KM=0
      INDEX9=0
      INDEX2=0
      GO TO 20
   21 KSTEP=4
      RETURN
    7 KSTEP=2
      RETURN
   20 KSTEP= 1
      RETURN
   15 KSTEP=3
      RETURN
      END
      FUNCTION GWCP(G,WCP)
      DIMENSION WCP(10,4)
      DO 2 I=1,10
      IF(WCP(I,1).LT.0.001) GOTO 3
   2  L=I
   3  IF(G.LT.WCP(1,1)) GOTO4
      DO 5 I=1,L
   5  IF(G.GT.WCP(I,1)) J=I
      IF(J.EQ.L) GOTO 6
      X=(WCP(J+1,2)-WCP(J,2))/(WCP(J+1,1)-WCP(J,1))*(G-WCP(J,1))+
     *WCP(J,2)
      GO TO 7
   4  X=(WCP(2,2)-WCP(1,2))/(WCP(2,1)-WCP(1,1))*(G-WCP(1,1))+
     *WCP(1,2)
      GO TO 7
   6  X=(WCP(L,2)-WCP(L-1,2))/(WCP(L,1)-WCP(L-1,1))*(G-WCP(L,1))+
     *WCP(L,2)
    7 GWCP=X
      RETURN
      END
      SUBROUTINE WAY1(M,INDEX2,GB,AN,PPI,SIGMA,DELTQ,DSIGMA,K,FX,FFX,PB)
      DIMENSION FX(24),FFX(24)
      COMMON/A23/INDEX6,INDEX7,INDEX8,INDEX9,GUIDE
      COMMON/A24/IL,IM,IN,KM,KMI
      COMMON/A25/SL1(7),SL2(6),DL1(11),DL2(10)
      COMMON/A25R/PDL1(11),PDL2(10)
      COMMON/A25V/SIGF
      COMMON/A26/PIM,PIKPI
      COMMON/A27R/SIGKKP,SGKP,SMKP
      COMMON/C2/KSI,SL21(6)/C2R/SL22(6)
      COMMON/A37/WAG1,WAG2,IND8,IND8A
      COMMON/A38/IND8R
      DIMENSION GB(10),AN(10,3)
 1111 FORMAT(20I6)
 2222 FORMAT(10F12.6)
      K=0
      IF(INDEX6.EQ.1) GO TO 20
      IF(IND8A.GT.0.AND.INDEX8.GE.1) INDEX8=5
      IF(IND8A.GT.0.AND.INDEX8.GE.1) DELTQ=ABS(DELTQ)
      IF(IND8.LT.100) GO TO 20
      AN(M,2)=AN(M,2)-DELTQ
      INDEX8=4
      GUIDE=KM+AN(M,2)
      K=1
      RETURN
   20 CONTINUE
      IF(INDEX2.EQ.0) GO TO 1
      IF(INDEX6.EQ.0) GO TO 2
      INDEX9=1
      IF(ABS(DSIGMA).LE.GB(10)) GO TO 13
      IF(DSIGMA.GT.0.) DSIGMA=-DSIGMA
      DSIGMA=+DSIGMA/2.
    3 K=1
      RETURN
   13 IF(DSIGMA.GT.0.) DSIGMA=-DSIGMA
      DSIGMA=DSIGMA
      GO TO 3
    2 IF(INDEX8.LT.4) GO TO 4
      DO 51 J=1,11
      DL1(J)=PDL1(J)
      IF(J.LE.10) DL2(J)=PDL2(J)
   51 CONTINUE
      SMKP=SGKP
      DO 60 J=1,24
   60 FX(J)=FFX(J)
      SIGMA=SIGF
      INDEX6=1
      INDEX7=0
      INDEX8=0
      DELTQ=-DELTQ
      GUIDE=10*IM+KM+SIGMA
      IM=0
      GO TO 5
    4 INDEX7=1
      MQ=0
    6 IF(MQ.EQ.0) GO TO 7
      IF(DELTQ.LT.0.) DELTQ=-DELTQ
      DELTQ=0.5*DELTQ
      GO TO 5
    7 IF(ABS(DELTQ).LE.GB(9)) GO TO 15
      IF(DELTQ.GT.0) DELTQ=-DELTQ
      DELTQ=+0.5*DELTQ
      GO TO 5
    1 IF(INDEX8.NE.1) GO TO 8
      IF(IND8R.EQ.0) DELTQ=-0.3*GB(8)
      IND8R=1
      PIM=PPI
    8 IF(INDEX6.EQ.0) GO TO 9
      PIKPI=PIM/PPI
      IF(INDEX9.EQ.0) GO TO 30
      IF(ABS(DSIGMA).GT.GB(10)) GOTO 31
      DO 50 J=1,11
      DL1(J)=PDL1(J)
      IF(J.LE.10) DL2(J)=PDL2(J)
   50 CONTINUE
      SMKP=SGKP
      DO 61 J=1,24
   61 FX(J)=FFX(J)
      SIGMA=SIGF
      INDEX9=0
      IF(SL22(6).LT.0.1) GO TO 100
      DO 101 LI=1,5
      SL21(LI)=SL22(LI)
      SL22(LI)=0
  101 CONTINUE
      SL22(6)=0.
      IF(KMI.EQ.1) KSI=1
  100 CONTINUE
      GUIDE=10*IM+KM+SIGMA
      IM=0
      GO TO 77
   31 CONTINUE
      IF(DSIGMA.LT.0.) DSIGMA=-DSIGMA
      DSIGMA=DSIGMA/2.
      GO TO 77
   30 CONTINUE
    9 IF(GUIDE.GE.5) GO TO 10
      IF(INDEX7.NE.1) GO TO 5
      IF(INDEX8.NE.0) GO TO 5
      MQ=1
      GO TO 6
    5 AN(M,2)=AN(M,2)+DELTQ
      IF(INDEX8.EQ.0) GO TO 10
      IF(INDEX8.EQ.5) GO TO 10
      INDEX8=INDEX8+1
      IF(INDEX8.EQ.4) IND8A=IM
      IF(INDEX8.NE.4) GO TO 10
      DELTQ=0.25*GB(9)
      AN(M,2)=GUIDE-KM+DELTQ
   10 IF(GUIDE.LE.5.) GO TO 77
      IF(INDEX9.EQ.1) GOTO 77
      IF(PIKPI.GE.PB) RETURN
   77 K=1
      RETURN
   15 IF(DELTQ.GT.0.) DELTQ=-DELTQ
      DELTQ=DELTQ
      GO TO 5
      RETURN
      END
      SUBROUTINE WAY2(STAGE,FX,GB)
      COMMON/A8/ALFA1/A9/RRR,AT1,UK2,UK4
      COMMON/A23/INDEX6,INDEX7,INDEX8,INDEX9,GUIDE
      COMMON/A24/IL,IM,IN,KM,KMI
      COMMON/A25/SL1(7),SL2(6),DL1(11),DL2(10)
      COMMON/A27/SIGMA,DSIGMA,DELTQ
      COMMON/A28/BET2,PPI
      COMMON/A29/AT2
      COMMON/A34/PM
      COMMON/A27R/SIGKKP,SGKP,SMKP
      COMMON/A36/IPEP,KKM,SGM(20),DHQMIN(20),INDEXN
      DIMENSION FX(24),STAGE(24)
      DIMENSION GB(10)
      KMI=GUIDE-10*IN
      IF(KMI.EQ.2) GO TO 2
      DO 3 J=1,7
    3 SL1(J)=DL1(J)
      BET2=DL1(8)
      GB(2)=DL1(9)
      GB(4)=DL1(10)
      PPI=DL1(11)
      GB(1)=PM
      SIGKKP=SMKP
      GO TO 4
    2 DO 5 J=1,6
    5 SL2(J)=DL2(J)
      GB(2)=DL2(7)
      AT2=DL2(8)
      GB(4)=DL2(9)
      PPI=DL2(10)
      GB(1)=PM
    4 DO 6 J=1,24
    6 STAGE(J)=FX(J)
      IF(IPEP.GT.0) RETURN
      SIGMA=SIGMA-DSIGMA
      GUIDE=GUIDE-DSIGMA
      RETURN
      END
      SUBROUTINE WAY3(INDEX3,IND1,IND2,IND3,IND4,HQ,DHQ)
      COMMON/A21/KRAN(40)
      IF(KRAN(32).EQ.1) WRITE(7,*) '   WAY3 START '
      IF(KRAN(32).EQ.1) WRITE(7,*) 'INDEX3,IND1,IND2,IND3,IND4,HQ,DHQ'
      IF(KRAN(32).EQ.1) WRITE(7,*) INDEX3,IND1,IND2,IND3,IND4,HQ,DHQ
      IF(INDEX3.NE.0) GO TO 53
   55 INDEX3=IND1+IND2
      GO TO 56
   53 IF(INDEX3.NE.1) GO TO 56
      IF(IND2.GT.0) GO TO 55
   56 IND3=IND1+IND2
      IF(INDEX3.NE.IND3) GO TO 64
      IF(IND4.EQ.1) GO TO 202
      IF(INDEX3.EQ.3) GO TO 62
      DHQ=0.005
   63 HQ=HQ-DHQ
      IF(KRAN(32).EQ.1) WRITE(7,*) 'INDEX3,IND1,IND2,IND3,IND4,HQ,DHQ'
      IF(KRAN(32).EQ.1) WRITE(7,*) INDEX3,IND1,IND2,IND3,IND4,HQ,DHQ
      IF(KRAN(32).EQ.1) WRITE(7,*) '   WAY3 FINISH'
      RETURN
   62 DHQ=-0.005
      GO TO 63
   64 IF(IND4.EQ.1) GO TO 201
      IND4=1
  200 DHQ=-DHQ/2.
      GO TO 63
  201 IF(DHQ.GT.0.) GO TO 200
   61 DHQ=DHQ/2.
      GO TO 63
  202 IF(DHQ.GT.0.) GO TO 61
      GO TO 200
      END
      ! 转子参数计算
      SUBROUTINE ROTOR(II,RK,HA,RKI,HAI,STAGE,HQ,QLA1,
     *C1U,C2A,C2U,C2,BETA2,ALFA2,LA2,AT2,ALFA4,LA4,C4A,
     *HZ,HZ1,HT,PIK,PIT,PI,AKPX,INDEX5)
      COMMON/A2/R1/A6/R2,DK1,DK2,DK4,DDK1,DDK2,DDK4
      COMMON/A7/UK1/A8/ALFA1
      COMMON/A9/RRR,AT1,UK2,UK4/A10/AKH
      COMMON/A11/C2AO,C4AO,ALFA4O/A13/AK
      COMMON/A16/SIGMAK,SIGMAA,DSITAK,DSITAA
      COMMON/A21/KRAN(40)
      DIMENSION RK(20),HA(20),RKI(10),HAI(10),STAGE(24)
      REAL LA2,LA4
  412 FORMAT('   HQ QLA1 C1A C1U ALFA1 BETA1 AM11 AL11
     *C2A C2U ALFA2 BETA2 LA2 C4A ALFA4 LA4
     *AT1 AT2 HZ HZ1 HT PIK PIT PI AKPD')
  500 FORMAT('   ROTOR   崁梹?')
  600 FORMAT('   ROTOR   妿崡厭')
  700 FORMAT(10F12.6)
      IF(KRAN(28).EQ.1) WRITE(7,500)
      C2AOA=HAI(9)
      C1AO=STAGE(1)
      AKP=RK(19)
      BTK=RK(7)
      HO=STAGE(5)
      AKPZ=STAGE(4)
      XFK=RK(11)
      AIOK=STAGE(16)
      B1O=RKI(4)
      B2O=RKI(2)
      A3O=HAI(4)
      A4O=HAI(2)
      AIOA=HAI(3)
      XFA=HA(11)
      HZO=STAGE(13)
      HZO1=HZO*9.81/(UK1*UK1)
      Y=(HQ-1.)/(HQ*HZO1/(R1*R1))
      YR=1.17-1.25*HZO1/(R1*R1)
      IF(YR-0.01)1,2,2
    1 YR=0.01
      GO TO 3
    2 IF(YR-0.7)3,3,4
    4 YR=0.7
    3 HQR=1./(1.-YR*HZO1/R1**2)
      IF(Y.LE.YR) RETURN
      INDEX5=1
      IF(YR-0.35)6,7,7
    6 ETA1=1.-0.5*YR**2
      GO TO 8
    7 ETA1=1.-0.5*YR**2-2.2*(YR-0.35)**2
    8 ETAK=1.
      ETAA=1.
      IF(XFK.GT.0.5) ETAK=20.*XFK-9.
      IF(XFA.GT.0.5) ETAA=20.*XFA-9.
      ETAC=(ETAK+ETAA)/2.
      ETA1=1.-(1.-ETA1)*ETAC
      RKPZ=AKPZ*ETA1
      CALL ENTRI(HQR,C1AO,UK1,ALFA1,AK,RRR,AT1,R1,
     *C1AR,QLA1R,B1R,AM1R,AL1R)   !?????????????
      AIK=RK(5)-B1R   !B1R入口气流角
      CALL HZH(HQR,AIK,AIOK,UK1,R1,BTK,
     *HO,AKH,AKP,AT1,AK,RRR,HT,HZ1,HZ,PIT)  !功与压比计算
      CALL PPZPKC(HT,RKPZ,QLA1R,C1AR,C2AOA,PIKR,ALFA2R,ALA2R,C2AR,C2UR,
     *AT2,PIR,PITR,ALFA4R,ALA4R,C4AR,II)
      SIGKR=PIKR/PITR
      SIGAR=PIR/PIKR
      CALL SIGM(B1O,B2O,AIK,AIOK,XFK,DSITAK,AL1R,AK,SIGK) 
      DELK=SIGKR-SIGK  !DSITAK转子总压损失系数
      AIA=HA(5)-ALFA2R
      CALL SIGM(A3O,A4O,AIA,AIOA,XFA,DSITAA,ALA2R,AK,SIGA)
      DELA=SIGAR-SIGA
      CALL ENTRI(HQ ,C1AO,UK1,ALFA1,AK,RRR,AT1,R1,
     *C1A,QLA1,B1,AM11,AL11)
      AIK=RK(5)-B1
      CALL SIGM(B1O,B2O,AIK,AIOK,XFK,DSITAK,AL11,AK,SIGK)
      SIGMRK=DELK+SIGK
      CALL HZH(HQ,AIK,AIOK,UK1,R1,BTK,
     *HO,AKH,AKP,AT1,AK,RRR,HT,HZ1,HZ,PIT)
      PIK=PIT*SIGMRK
      CALL BEXIT(HZ,HZ1,C1A,PIK,QLA1,  !!!!!!!!!!!!!!
     *C1U,C2A,C2U,C2,BETA1,BETA2,ALFA2,LA2,AT2)  
      AIA=HA(5)-ALFA2
      CALL SIGM(A3O,A4O,AIA,AIOA,XFA,DSITAA,LA2,AK,SIGA)
      SIGMRA=SIGA+DELA
      PI=PIK*SIGMRA
      ! 计算出口流量及流速
      CALL AEXIT(C2A,AT2,PI,QLA1,C2AOA,ALFA4,LA4,C4A,II)
      SMUL=1.
      IF(HZ.LT.0.) SMUL=-1.
      IF(ABS(AT1-AT2).LT.1.E-5) AT2=AT1*(1.+SMUL*1.E-5)
      AKPX=(PI**((AK-1.)/AK)-1.)/(AT2-AT1)*AT1
      IF(KRAN(28).EQ.0) RETURN
      WRITE(7,412)
      WRITE(7,700) HQ,QLA1,C1A,C1U,ALFA1,BETA1,AM11,AL11,
     *C2A,C2U,ALFA2,BETA2,LA2,C4A,ALFA4,LA4,
     *AT1,AT2,HZ,HZ1,HT,PIK,PIT,PI,AKPX
      WRITE(7,600)
      RETURN
      END
      SUBROUTINE SIGM(B1O,B2O,AIK,AIOK,XFK,DSITAK,AL1R,AK,SIGK)
      COMMON/A21/KRAN(40)
      REAL AL1R
  411 FORMAT('   DSI AKS SIGK')
  501 FORMAT('   SIGM   崁梹?')
  601 FORMAT('   SIGM   妿崡厤')
  700 FORMAT(6F20.6)
      IF(KRAN(25).EQ.1) WRITE(6,501)
      DBO=B2O-B1O
      DI=AIK-AIOK
      DIR=ABS(DI/DBO)
      IF(DIR-2.)1,1,2
    1 DSI=1.+(3.*DIR+5.)*DIR
      GO TO 3
    2 DSI=16.*DIR-9.
    3 IF(XFK-0.5)4,4,5
    4 AKS=1.
      GO TO 6
    5 AKS=20.*XFK-9.
    6 DSIR=DSITAK*(1.+(DSI-1.)*AKS)
      FF=Q(AL1R,AK)*AL1R*(2./(AK+1.))**(1./(AK-1.))
      SIGK=1.-AK/(AK+1.)*FF*DSIR
      IF(KRAN(25).EQ.0) RETURN
      WRITE(7,411)
      WRITE(7,700) DSI,AKS,SIGK
      WRITE(7,601)
      RETURN
      END
      SUBROUTINE HZH(HQR,AIK,AIOK,UK1,R1,BTK,HO,AKH,AKP,
     *AT1,AK,RRR,HT,HZ1,HZ,PIT)    !主特性理论能头计算
      COMMON/A21/KRAN(40)
      COMMON/A31/F1,XFK,XFA
  403 FORMAT('   H HOB DI1R A HZ HT PIT')
  503 FORMAT('   HZH   崁梹?')
  603 FORMAT('   HZH   妿崡厤')
  700 FORMAT(6F20.6)
      IF(KRAN(27).EQ.1) WRITE(7,503)
      HXFK=1.
      IF(XFK.GT.0.5) HXFK=1.5*XFK+0.25
      IF(HQR-1.3)1,1,2
    1 H=0.9*(1-HQR)
      GO TO 8
    2 H=1.03-HQR
    8 HOB=H*AKP   !AKP压头特性校正系数
      DI1R=(AIK-AIOK)/SQRT(BTK)
      DI1R=DI1R*57.296
      IF(DI1R+10.)3,3,4
    3 IF(DI1R+50)5,6,6
    5 A=1.28
      GO TO 7
    6 A=0.93-0.007*DI1R
      GO TO 7
    4 A=1.
    7 H=HOB*A
      H=H*HXFK
      HT=H*R1**2+HQR*HO
      HZ1=HT*AKH
      HZ=HZ1*UK1**2/9.81
      PIT=(1.+HZ/(AK*RRR*AT1/(AK-1.)))**(AK/(AK-1.))   !压比
      IF(KRAN(27).EQ.0) RETURN
      WRITE(7,403)
      WRITE(7,700) H,HOB,DI1R,A,HZ1,HZ,HT,PIT
      WRITE(7,603)
      RETURN
      END
      SUBROUTINE ENTRI(HQR,C1AO,UK1,ALFA1,AK,RRR,AT1,R1,
     *C1AR,QLA1R,B1R,AM1R,AL1R)
      COMMON/A21/KRAN(40)
      REAL LA1R
  402 FORMAT('   C1AR QLA1R B1R AM1R AL1R')
  502 FORMAT('   ENTRY   崁梹?')
  602 FORMAT('   ENTRY   妿崡厤')
  700 FORMAT(6F20.6)
      IF(KRAN(26).EQ.1) WRITE(7,502)
      C1AR=HQR*C1AO
      LA1R=C1AR*UK1/(SIN(ALFA1)*SQRT(2.*9.81*AK*RRR*AT1/(AK+1.))) 
      QLA1R=Q(LA1R,AK)
      CALL AL1(C1AR,R1,B1R,AM1R,AL1R)  !AL1R相对速度数
      IF(KRAN(26).EQ.0) RETURN
      WRITE(7,402)
      WRITE(7,700) C1AR,QLA1R,B1R,AM1R,AL1R
      WRITE(7,602)
      RETURN
      END
       FUNCTION COD(A,N,B,M,K)
        DOUBLE PRECISION COD,C,D,A3
        INTEGER*4 I1,I2
        COMMON/A21/KRAN(40)
  999  FORMAT(F15.8,1I7,F15.8,2I7)
      IF(KRAN(1).EQ.1) WRITE(15, 999) A,N,B,M,K
1000  FORMAT(5I7,2F15.8,3F20.1)
      A1=A*10.**(N+1)
      I1=IFIX(A1)   !向下取整
      A2=A*10.**N
      I2=IFIX(A2)
      I2=I2*10
      I3=I1-I2
      J=0
      IF(I3.GE.5) J=10
      I1=(I2+J)/10
      A3=DFLOAT(I1)    !将I1转换为双精度实数
      C=A3*10**K
      B1=B*10.0**(M+1)
      I1=IFIX(B1)
      B2=B*10.0**M
       I2=IFIX(B2)
       I2=I2*10
       I3=I1-I2
       J=0
      IF(I3.GE.5) J=10
      I1=I2+J
      I2=I1/10
       D=DFLOAT(I2)
      COD=C+D
      RETURN
      END
       SUBROUTINE HARCK(NZAP,GR,PIR,KPDR,NR,LA3,TB,
     *                 K,D2,BT2,K1,K2,K3,R,U2,
     *                 KK,GI,PII,KPDI,NOM)
      REAL KPDR,K,LA3,N(10),KPD(10),NO(10),KPDI(20),KPDO,KPOP(10),
     *KPOL(10),KPO(20),NZAP,NR
      DIMENSION G(10),U2(10),GI(20),HTI(20),VP(10),VL(10),V(20),
     *HP(10),H(20),HT(10),PII(20),NTO(20)
      IGG=0
      KB=1
      KDOP=0
      KB1=KB+KDOP
      N(KB)=NZAP
    1 FORMAT (6I2)
    2 FORMAT (8F8.4)
    3 FORMAT (10F7.1)
   17 FORMAT ('0',20X,'垜晭剭泤  剙崓泤'/' ',20X,16('_')/
     *'0',3X,'梹憭帓€ 倫€檯崍?N, 巵/寛?',14('.'),F9.3/
     *' ',3X,'彁垈厔厤崨?悁憰巹 値噭摃€ G 彁, 妰/C ',F9.3/
     *' ',3X,'憭厪厤?弾倹槄崍?剙€倠厤垷 ',14('.'),F9.3/
     *' ',3X,'妿潝攬枅厤?弾媴噸巸?剠墤拏垷 ',10('.'),F9.3/
     *' ',3X,'彁垈厔厤崁?憡帎帒挏 弲悈?剤敂搰帎帉 ',F9.3/
     *' ',3X,'拝審厫€挀悁 崁 倳巹? TO, K ',13('.'),F9.3/
     *' ',3X,'妿潝攬枅厤?垏帩崚悗張梾憡巸?彁帠厬憖 K ',F9.3/
     *' ',3X,'剤€寘拹 妿媴憖 崁 倹晭剠  D2, MM ',8('.'),F9.3/
     *' ',3X,'搩帇 嫀弨拪?崁 倹晭剠  B2L, 儛€創?',5('.'),F9.3/
     *' ',3X,'崕寘悁 姁垈洉: K1=',I1,'  K2=',I1,'  K3=',I1)
      R=287.023
      HAR=K/(K-1)*R*TB*(PIR**((K-1)/K)-1)
      HTR=HAR/KPDR
      DO 4 I=1,KB1
      U2(I)=3.1416*D2*N(I)/60.
    4 CONTINUE
      UR=3.1416*D2*NR/60.
      G(KB)=GR
      KPD(KB)=KPDR
      HT(KB)=HTR/UR**2
    5 FORMAT (' ',3X,50('*'))
      NTO(1)=1
      DO 7 I=2,20
    7 NTO(I)=NTO(I-1)+1
      DO 8 I=1,KB1
      NO(I)=N(I)/NR
      ON=NO(I)
      CALL REGIM (ON,K1,K2,K3,BT2,PIR,K,GO,HTO,KPDO)
      G(I)=GO*G(KB)
      HT(I)=HTO*HT(KB)
      KPD(I)=KPDO*KPD(KB)
    6 CONTINUE
  300 FORMAT(//10X,'帓崕憟拝嫓崁?梹憭帓€ 倫€檯崍? N=',F6.2)
   18 FORMAT ('0',11X,'PEVIM',I1,'   梹憭帓€ 倫€檯崍?  N=',F9.3/' ',
     *20X,'帄悡啀€?憡帎帒挏    U2=',F7.3/' ',10X,
     *40('-')/' ',12X,'拵棅? !  G PP   !    P    !   KPD   !'/' ',
     *10X,40('-'))
      CALL PRAVO (ON,LA3,BT2,KPOP,VP,HP,VMA)
      GMA=VMA*G(I)
      GMI=DPZ(LA3,ON,GMA)
      VMI=GMI/G(I)
      KL=5
      KP=11
      KK=KP+3
      KR=KP+1
      DO 14 J=6,KP
      V(J)=VP(J-KL)
      H(J)=HP(J-KL)
      KPO(J)=KPOP(J-KL)
      GI(J)=V(J)*G(I)
      HTI(J)=H(J)+V(J)*HT(I)
      KPDI(J)=KPO(J)*KPD(I)
      HAI=KPDI(J)*HTI(J)*U2(I)**2
      PII(J)=(HAI*(K-1.)/K/R/TB+1.)**(K/(K-1.))
   14 CONTINUE
      CALL LEVO1(ON,VMI,LA3,VL,KPOL)
      DO 11 J=1,KL
      V(J)=VL(J)
      KPO(J)=KPOL(J)
      GI(J)=V(J)*G(I)
      KPDI(J)=KPO(J)*KPD(I)
      PII(J)=-0.1*ON*(GI(J)-GI(6))+PII(6)
   11 CONTINUE
      PII(KK)=1.
      DPI=(PII(KP)-1.0)/3
      PII(KP+1)=PII(KP)-DPI
      PII(KP+2)=PII(KP)-2*DPI
      DO 20 J=KR,KK
      HAI=K/(K-1)*R*TB*(PII(J)**((K-1)/K)-1)
      HTI(J)=H(KP)+V(KP)*HT(I)
      KPDI(J)=HAI/HTI(J)/U2(I)**2
   20 GI(J)=GI(KP)
   19 FORMAT (' ',13X,I2,4X,F9.3,1X,F9.3,1X,F9.3)
      IF(IGG.EQ.0) GOTO 8
  400 FORMAT(4F10.4)
    8 CONTINUE
 200  RETURN
      END
      SUBROUTINE REGIM (N,K1,K2,K3,BL,PI,K,G,H,KP)
      REAL N,KP,K,NB,N1
      N1=N
      GO TO (1,2), K1
    1 A=0.5867*(PI**((K-1)/K)-1)+0.9933
      G=(-A+SQRT(A**2*(1-4*N**2+4*N)+4*A*N*(2*N-3)+4*N*(2-N)))/2/(1-A)
      GOTO 3
    2 RESERV=0
    3 CONTINUE
      GO TO (4,5), K2
    4 H=(-0.3125)*N**2+0.375*N+0.9375
      GOTO 7
    5 REZERV=0
    7 CONTINUE
      GO TO (14,6,11,12),K3
    6 IF(N.LT.0.9) N=0.9
   14 BE=BL/57.3
      IF(N.LT.1.0) GOTO 16
      Z=0.2101*BE**2-0.2933*BE+0.1024
      A=-(0.02)/Z
      B=-(0.0183*BE-0.0528)/Z
      C=(0.2101*BE**2-0.275*BE+0.0696)/Z
      KP=A*N**2+B*N+C
      GOTO 13
   11 KP=-0.1224*N**2+0.0979*N+1.03042
      GOTO 13
   16 CONTINUE
      BN=BE*57.3
      BN=(BN-70.)/10.
      NB=( N-0.7)/0.15
      KP=1.017179-.003982*BN-.001196*NB+.000438*BN**2+
     *.002203*BN*NB-.00408*NB**2
   12 REZERV=0
   13 CONTINUE
      N=N1
      RETURN
      END
      SUBROUTINE PRAVO (ON,LA,B,KP,V,H,VM)
      REAL LA,KP(10)
      DIMENSION V(10),H(10)
      K=5
      KT=K+1
      SL=LA
      R=(0.0461*B/57.3+0.0096)*LA**(-1.112)
      IF(ABS(ON-1.).LE.0.0001) GOTO 2
      RO=EXP(2.773*(1-ON))
      R=R*RO
      LA=(-0.3857)*ALOG(R)+0.2243
    2 VM=R+1-0.0005
      V(1)=1.
      H(1)=0.
      KP(1)=1.
      S=(VM-V(1))/K
      DO 1 I=2,KT
      V(I)=V(I-1)+S
      H(I)=-0.0001*V(I)**(58.4*(2*LA-1))-0.78125*(V(I)-1)
      KP(I)=1-R*(1+0.5476*(1-(V(I)-1)/R)-(2.395-1.6949*(V(I)-1)/R-
     *0.7002*((V(I)-1)/R)**2)**0.5)
    1 CONTINUE
      LA=SL
      RETURN
      END
      SUBROUTINE LEVO1(ON,VMI,LA,V,KP)
      REAL LA,KP(10)
      DIMENSION V(10)
      K=5
      KT=K+1
      V(KT)=1.
      I=KT
      S=(V(KT)-VMI)/K
      SL=LA
      OL=(0.725*ON+0.275)
      LA=OL*SL
    1 IF (I.EQ.1) GOTO 2
      V(I-1)=V(I)-S
      KP(I-1)=1.-(ABS(2.8333*LA-0.9833)*(1.-V(I-1))**2+(0.2333*LA-
     *0.1133)*(1.-V(I-1)))
   4  I=I-1
      GO TO 1
    2 CONTINUE
      LA=SL
      RETURN
      END
      FUNCTION DPZ (LA,ON,G)
      REAL LA
      P=(-0.8307)*LA**2+1.7971*LA-0.1115
      IF (ON.GE.1.) GOTO 1
      PO=1.-0.8672*(1.-ON)**1.2188
      P=P*PO
    1 DPZ=P*G
      RETURN
      END
      SUBROUTINE OPTIMA(II,RK,HA,RKI,HAI,
     *GB,GKA,X,Y,GDK,GDA,STAGE,INDEX1,INDEX4,OPTIMS,GKI,M)
      DIMENSION OPTIMS(10),GKI(10)
      DIMENSION RK(20),HA(20),RKI(10),HAI(10),GB(10),GKA(4,30),
     *X(16),Y(16),GDK(6,30),GDA(6,30),GD(6),GK(4),STAGE(24)
      COMMON/A22/C(4)
      COMMON/A20/ARKI(30,10),AHAI(30,11)
      COMMON/A7/UK1/A8/ALFA1/A9/RRR,AT1,UK2,UK4/A10/AKH
      COMMON/A16/SIGMAK,SIGMAA,DSITAK,DSITAA
      COMMON/BLOK/AP1,A1,R,UKR
      DO 1414 JJ=1,10
      HAI(JJ)=0.
1414  RKI(JJ)=0.
      DO 5 K=1,4
    5 GK(K)=GKA(K,II)
      IF(INDEX4-1)7,6,6
    7  CALL GEOM (RK,X,Y,GD,HA(19))   ! 计算动叶栅几何参数
      DO 3 K=1,6
    3 GDK(K,II)=GD(K)
      CALL GEOM(HA,X,Y,GD,HA(20))     ! 计算静叶栅几何参数
      DO 4 K=1,6
    4 GDA(K,II)=GD(K)
    6  UK2=UK1*RK(2)/RK(1)            ! 动叶进口速度
      UK4=UK1*HA(2)/RK(1)             ! 静叶进口速度
      CALL POINT(RK,HA,GK,GB,GDK,
     *GDA,HAI,RKI,STAGE,II,INDEX1,OPTIMS,M)
      IF(ABS(GK(1)-GKA(1,II)).LT.0.01) GO TO 70
      IF(GK(1).LT.GKI(II)) GKI(II)=GK(1)
   70 CONTINUE
   71  FORMAT(8F9.4)
      RETURN
      END
      ! 最佳入口攻角计算 PDF - P17,2.入口相对气流角及最佳攻角的确定
      SUBROUTINE AI0(RK,DI,BETA2,AI0K,BETA1,FER)
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
      ! 转子效率计算
      FUNCTION AKPD(RK,AL10)
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
      ! 根据流道尺寸对效率的修正值（P20,3-23）
      FUNCTION HP(AHP,AL10)
      HR=AL10**2*(25.-AHP)/144.
      IF(HR.LT.0) HR=0.
      HP=HR
      RETURN
      END
      ! 落后角及出口相对气流角的确定
      SUBROUTINE DELTAK(RK,GK,J,E,C2C1,B20,DELT,BETA2)
      DIMENSION RK(20),GK(4)
      P=57.296
      B20=B20*P
      RK(6)=RK(6)*P
      E=E*P
      IF(J-1)1,1,2
    2 IF(RK(16))10,11,10
   10 XX=1.
    6 XC=0.455/(C2C1*RK(2)/RK(1)-0.545)
      GOTO7
   11 XX=1.1
      IF(RK(3)-0.4)3,3,4    !RK(3)相对轮毂直径，即轮毂比
    3 XC=1.
      GOTO7
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
      ! 效率修正系数Ok计算
      SUBROUTINE OKPD(RK,C1A0,C1A0K,AL10,ANW,OK,FK)
      DIMENSION RK(20)
      AL=AL10        !进口无因次速度数
      ANW=C1A0/C1A0K
      IF(RK(16))1,2,1
    2 IF(AL10-0.5)3,3,4
    3 AL10=0.5
    4 FK=AL10**2.-0.97*AL10+0.3
      GOTO7
    1 IF(AL10-0.9)5,5,6
    5 AL10=0.9
    6 FK=0.46*AL10-0.35
    7 OK=1-10.*(ABS(ANW-1.)**1.334)*FK  !转子效率修正系数
      AL10=AL
      RETURN
      END
      SUBROUTINE GEOM (PK,X,Y,GD,GAGT)        !确定叶轮几何参数
      DIMENSION PK(20),X(16),Y(16),GD(6)
      YMAX=0.
      DO 18 L=1,16
        IF(Y(L).GT.YMAX) YMAX=Y(L)
   18 CONTINUE
      DK1=PK(1)           ! 叶栅进口外缘直径
      DK2=PK(2)           ! 叶栅出口外缘直径
      D1O=PK(3)           ! 叶栅进口处轮毂比
      D2O=PK(4)           ! 叶栅出口处轮毂比
      B1K=PK(5)           ! 叶栅几何进气角B1k
      B2K=PK(6)           ! 叶栅几何出气角B2k
      B1T=PK(7)           ! 叶栅稠度
      B=PK(8)             ! 叶片弦长
      L=PK(10)            ! 亚音速或超音速
      XF=PK(11)
      C=PK(12)
      AG2=PK(15)
      G12=PK(16)
      B1SP=B1K-G12
      E=B2K-B1K
      IF(AG2.NE.0) GO TO 50
      A=XF/(1.-XF)
      QE=0.5*(1./A+A)/TAN(E)   !?????
      HI1=ATAN(1./(A*(QE+SQRT(1.+QE**2)))) !前缘角
      HI2=E-HI1   !后缘角
      V=E/2.
      KP=0        ! KP为0表示最大弯度在叶形中间
      IF(ABS(XF-0.5).GT.0.0001) KP=1
      T=B1K+HI1    !安装角
      IF(Y(1).GT.0..OR.Y(16).GT.0.) GOTO 5   !叶栅喉部尺寸
      Y(1)=Y(2)-X(2)*(Y(3)-Y(2))/(X(3)-X(2))
      Y(16)=Y(15)-(X(15)-1.)*(Y(14)-Y(15))/(X(14)-X(15))
    5 NT=30
      S=1./(NT-1)
      A=TAN(HI1)
      B=TAN(HI2)
      R=(.5/SIN(V))**2
      AG2=10.
      DO 15 J=1,NT
      X1=(J-1)*S
      CALL PARINV(X1,X,Y,16,YT)    !插值  
      Y1=YT*C/YMAX*.5
      IF(KP-1)6,7,8
    6 H=X1*(1.-X1)
      HI=(.5-X1)/SQRT(.25/A**2+H)
      YC1=SQRT(R-.25+H)-SQRT(R-.25)
      GO TO 8
    7 H=1.-X1
      HI=A*B*(H**2*B-X1**2*A)/(X1*A+H*B)**2
      YC1=X1*H*A*B/(H*B+X1*A)
    8 YHI1=ATAN(HI)
      AG1=10.0
      DO 14 I=1,NT
      X2=(I-1)*S
      CALL PARINV(X2,X,Y,16,YT)
      Y2=YT*C/YMAX*.5
      IF(KP-1)9,10,11
    9 H=X2*(1.-X2)
      HI=(.5-X2)/SQRT(.25/A**2+H)
      YC2=SQRT(R-.25+H)-SQRT(R-.25)
      GO TO 11
   10 H=1.-X2
      HI=A*B*(H**2*B-X2**2*A)/(X2*A+H*B)**2
      YC2=X2*H*A*B/(H*B+X2*A)
   11 YHI2=ATAN(HI)
      AG=B1T*(SQRT((X1+COS(T)/B1T-X2)**2+(YC1+SIN(T)/B1T-YC2)**2)
     *-Y1*COS(YHI1)-Y2*COS(YHI2))     ! 喉部相对面积
      IF(AG.LT.AG1) GO TO 14
      AGT=AG1          !叶栅喉道尺寸与步长之比
      GO TO 40
   14 AG1=AG
   40 IF(AGT.GT.AG2) GO TO 50
   15 AG2=AGT
   50 R1=SQRT((1.+D1O**2)/2.)
      AGT=AG2
      IF(GAGT.GT..001)AGT=AGT*GAGT    ! 喉部面积
      R2=SQRT((1.+D2O**2)/2.)
      F2F1=((1.-D2O**2)*DK2**2)/((1.-D1O**2)*DK1**2)
      GD(1)=E     ! 气流折转角
      GD(2)=AGT   ! 喉部面积
      GD(3)=R1    ! 
      GD(4)=R2
      GD(5)=F2F1  ! 叶栅进出口面积比
      GD(6)=B1SP  ! 进气角与叶型楔角一半的差（可能用于计算攻角损失）
      RETURN
      END
      ! 动力粘性系数系数
      FUNCTION AMQ(T1)
      AMQ=1.7216*1E-6*(T1/273.)**1.5*395./(T1+122.)
      RETURN
      END
      ! 轴向相对速度修正系数计算
      FUNCTION ONWK(RK,ANWK,ANWA,TAW,AKPDK,AKPDA,FK,FA)
      DIMENSION RK(20)
      X=FA/FK
      IF(RK(16))1,3,1
    1 X=0.7*X
    3 AO=AKPDA*X*(1.-TAW)/(AKPDK*TAW)
      ONWK=(1.+AO*ANWA/ANWK)/(1.+AO*(ANWA/ANWK)**2.)
      RETURN
      END
      ! 攻角修正系数计算
      FUNCTION ATAKA(RK,AI0K,R1,C1A0K,B1,AL10)
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
      COMMON/A34/PM
      COMMON/A35/IJ
      COMMON/A36/IPEP,KKM,SGM(20),DHQMIN(20),INDEXN
      COMMON/A37/WAG1,WAG2,IND8,IND8A
      COMMON/A38/IND8R
      COMMON/A40/LA4
      COMMON/C1/III,INDEXC,HQKP,HQKPA
      COMMON/C3/HQK,HQA
      COMMON/STUP/IQPN
C      common/add/jicanshu
        DOUBLE PRECISION STUP,COD
  112 FORMAT('N STUP',4X,4HQLA1,8X,2HPI,7X,3HKPD,8X,2HHT,8X,2HT2,
     17X,3HC1A,8X,2HHQ,6X,4HHQCK,7X,3HHQC,5X,5HHQKPK,4X,5HHQKPA)
  113 FORMAT(2X,I2,1X,11F10.4)
  933 FORMAT(6F20.6)
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
      R1=SQRT((1+D1**2)/2)    ! ?????
      R2=SQRT((1+D2**2)/2)    ! 
      HRK=RK(9)               ! 转子叶片展弦比
      XFK=RK(11)              ! 转子叶片中线的最大挠度相对坐标
      XFA=HA(11)              ! 静子叶片中线的最大挠度相对坐标
      CK=RK(12)               ! 转子叶片型面最大相对厚度
      CA=HA(12)               ! 静子叶片型面最大相对厚度
      YY=RK(8)                ! 工作轮叶片弦长
      H1=RK(9)                ! 工作轮展弦比
      HB3=HA(9)               ! 导向器展弦比
      GG1=RK(7)   !工作轮叶栅稠度
      GG3=GG1*R1*DK1/(R2*DK2)  !?????????????????
      H2=DK2*(1-D2)/2         ! 转子出口流道径向高度的一半
      T1=YY/GG1   !动叶进口节距
      T2=YY/GG3   !动叶出口节距??????????????
      HB4=HB3*DK4*(1-DDK4)/(DK2*(1-DDK2))   !???????????
      BT3=HA(7)     !静叶叶栅稠度
      BT4=BT3*SQRT((1.+DDK2**2.)/
     *(1.+DDK4**2.))*DK2/DK4
      IF(RK(16))4449,4448,4449
 4448 UPR2=0.1
      GO TO 4450
 4449 UPR2=2
 4450 CONTINUE
      AT1=GB(2)                           ! 气体总温
      RRR=GB(6)                           ! 气体常数R
      ALFA1=GB(4)                         ! 工作轮进口处气流角
      ALFA1=ALFA1/PA                      ! 换算成弧度制
      A1=ALFA1
      F1=3.14159*DK1**2.*(1.-D1**2.)/4.   ! 
      AM=SQRT(AK*(2/(AK+1))**((AK+1)/(AK-1)))      !求流量系数K=AM/SQRT(R)
      IF(II-IL)224,224,225
  224 AA=AM*F1/SQRT(AT1)*SIN(ALFA1)
  225 AA=AA
      UK2=UK1*DK2/DK1
      UK4=UK1*DK4/DK1
      QLA=AN(M,2)         ! 流量系数
      ANT=AN(M,1)         ! 转速百分比
      AKG=GKA(2,II)       ! 空气流量储存系数
      AKH=GKA(3,II)       ! 理论压头减少系数
      IF(KRAN(5).EQ.0) GO TO 66
      WRITE(7,933) DK1,DK2,DK4,D1,D2,DDK1,DDK2,
     *DDK4,R1,R2,HRK,XFK,XFA,CA,CK,YY,T1,T2,
     *H1,HB3,HB4,BT3,BT4,UPR2,AT1,ALFA1,RRR,AK,AKP,AKH,AKG,UK1,UK2,UK4
      WRITE(7,933) QLA,ANT,UK1
   66 C1AO=STAGE(1)
      HTO=STAGE(5)
      AML1=STAGE(6)
      C2AO=STAGE(9)
      ADST=STAGE(4)
      B2O=STAGE(15)
      C4AO=STAGE(19)
      ALFA4O=STAGE(2)
      RE=STAGE(17)        ! 雷诺数?
      PPI=STAGE(22)
      QLA1=QLA*AKG*AA*SQRT(AT1)/(AM*F1*PPI*SIN(ALFA1))    !QLA-压气机入口折合流量
      AKR=SQRT(2*9.81*RRR*AT1*AK/(AK+1))                  !压气机入口处临界声速
      AL1=RLQMDA(QLA1)                                    !入口处无因次速度数
      C1A=AL1*AKR*SIN(ALFA1)/UK1                          !入口处速度系数
      HQ=C1A/C1AO
      CALL CRITIC(II,HQ,QLA1,C1AO,B2O,HTO,ADST,KAG,HAG,RE,UPR2,M,
     *RK,HA,RKI,HAI,STAGE,AN,GB,GKA,
     *INDEX2,INDEX3,C1U,C2A,BETA1,BETA2,ALFA2,LA2,AT2,
     *LA4,C4A,PI,PIK,PIT,HT,AKPX)         ! 临界？
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
      DO 211 J=1,10
      DL2(J)=0
  211 DL1(J)=0
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
      GB(2)=AT2
      GB(4)=ALFA4
      PPI=PPI*PI
      GB(1)= PI*GB(1)
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
C      IF(jicanshu.EQ.0) RETURN
	RETURN
      BETA2=BETA(C2A,R2,ALFA2)
      WRITE(7,113)II,QLA1,PI,AKPX,(beta1-rk(5))*180/3.14,
     *(beta1-rk(5))*180/3.14,(ha(6)-alfa4)*180/3.14
      RETURN
      END
      ! 分离参数计算
      SUBROUTINE BLOCK(II,HQ,QLA1,C1AO,B2O,HTO,ADST,RE,UPR2,
     *RK,HA,RKI,HAI,STAGE,
     *ALFA2,C2A,AT2,PI,PIK,PIT,LA4,C4A,HT,AKPX)
      DIMENSION RK(20),HA(20),RKI(10),HAI(10),STAGE(24)
      REAL LA2,LA4
      COMMON/A1/D1/A2/R1/A3/YY,H1,H2,T1,T2,CK
      COMMON/A5/HRK
      COMMON/A6/R2,DK1,DK2,DK4,DDK1,DDK2,DDK4
      COMMON/A7/UK1/A8/ALFA1/A9/RRR,AT1,UK2,UK4/A10/AKH
      COMMON/A11/C2AO,C4AO,ALFA4O/A12/AKP/A13/AK
      COMMON/A14/HB3,HB4,CA,BT3,BT4
      COMMON/A15/AKR/A43/FGK,FGKA
      COMMON/A16/SIGMAK,SIGMAA,DSITAK,DSITAA
      COMMON/A17/QLA4,ALFA4,QQLA1
      COMMON/A21/KRAN(40)
      COMMON/A22/C(4)/R/INDEX5
      COMMON/A28/BET2,PPI
      COMMON/A30/C1A,AML1
      COMMON/A31/F1,XFK,XFA/A32/ANT,AKG
      COMMON/A33/QLA1C,HQCK,HQC
  933 FORMAT(6F20.6)
 9992 FORMAT(3X,15H分离的级2X,2HN=,I2,4HFGK=,F8.4,2X,
     *5HFGKA=,F8.4)
      C2AOA=HAI(9)
      INDEX5=0
      AKP=RK(19)   !压头特性校正系数
      CALL YPPKC(HQ,C1AO,HTO,B2O,AML1,ADST,PIK,ALFA2,LA2,C2A,C2U,AT2,
     *II,RK,HA,RKI,HAI,STAGE,
     *AKP,XFK,XFA,UPR2,PI,PIT,ALFA4,LA4,C4A,HT,AKPX)
      IF(HQC.GT.0.) GO TO 2114
      IF(RK(20).GT.1E-6) GO TO 1000   ! 输入参数中提供：RK(20)分离边界上的空气动力载荷判据值
      CALL FGCA(ANT,RE,FGK,FGKA)      ! 计算分离判据
      GO TO 1001
 1000 FGK=RK(20)
      FGKA=FGK*1.2
 1001 IF(INDEX5.GT.0) RETURN
       XXX=AK*RRR*9.81
      HQCK=HQCP(C1AO,B2O,XXX,FGK)
      CALL OGCC(HQ,HQCK,C1AO,HTO,B2O,AML1,ADST,FGKA,
     *II,RK,HA,RKI,HAI,STAGE,
     *AKP,XFK,XFA,UPR2,HQC,ALIN)
 2114 IF(HQC-50.)2112,2112,2113       ! 到2113为计算分离区参数
 2112 IF(KRAN(39).EQ.1) WRITE(7,9992)II,FGK,FGKA
      IF(KRAN(39).EQ.1) WRITE(15,9992)II,FGK,FGKA
      CALL PHOY(HQC,UPR2,HTO,C1AO,B2O,H)
      HC=H*PHMK(AKP,XFK)
      CALL PHCY(HC,HQ,HQC,H)
      HT=H*R1**2+HQ*HTO
      CALL PKPDOY(HQC,UPR2,AML1,HAI(10),VTAC)
      CALL PKPDK(XFK,XFA,HQC,VTAC)
      HQT=HQ/HQC
      CALL PKPDCY(HQT,VTAT)
      VTA=VTAC*VTAT
      AKPX=VTA*ADST
      C1AC=HQC*C1AO
      QLA1C=QQL(C1AC)
      IF(KRAN(6).EQ.0) GO TO 2117
      WRITE(7,933) HC,H,HT,VTAC,VTAT,AKPX
 2117 CALL PPZPKC(HT,AKPX,QLA1,C1A,C2AOA,PIK,ALFA2,LA2,C2A,C2U,AT2,
     *PI,PIT,ALFA4,LA4,C4A,II)
      IF(KRAN(6).EQ.0) GO TO 2113
      WRITE(7,933) HQ,C1A,PIK,ALFA2,LA2,C2A,
     *C2U,AT2,PI,PIT,ALFA4,LA4,C4A,HT,AKPX
 2113 CONTINUE
      RETURN
      END
      SUBROUTINE OGCC(HQ,HQCK,C1AO,HTO,B2O,AML1,ADST,FGKA,
     *II,RK,HA,RKI,HAI,STAGE,
     *AKP,XFK,XFA,UPR2,HQC,ALIN)
      REAL LA2,LA4
      COMMON/A2/R1/A6/R2,DK1,DK2,DK4,DDK1,DDK2,DDK4
      COMMON/A7/UK1/A8/ALFA1
      COMMON/A9/RRR,AT1,UK2,UK4/A10/AKH
      COMMON/A11/C2AO,C4AO,ALFA4O/A13/AK
      COMMON/A16/SIGMAK,SIGMAA,DSITAK,DSITAA
      COMMON/A21/KRAN(40)
      DIMENSION RK(20),HA(20),RKI(10),HAI(10),STAGE(24)
  571 FORMAT('O G C C 崁梹拵')
  572 FORMAT('O G C C 妿崡厤?')
 1357 FORMAT(6F20.6)
      IF(KRAN(22).EQ.1) WRITE(7,571)
      A=1.
      B=0.
      HQT=HQ
      ALIN=1
      IF(HQT-HQCK)711,712,712
  711 CALL YPPKC(HQCK,C1AO,HTO,B2O,AML1,ADST,
     *PIK,ALFA2,LA2,C2A,C2U,AT2,
     *II,RK,HA,RKI,HAI,STAGE,
     *AKP,XFK,XFA,UPR2,PI,PIT,ALFA4,LA4,C4A,HT,AKPX)
      FGZ=FGZA(ALFA2,ALFA4,LA2)
      IF(FGZ-FGKA)714,714,715
  714 HQC=HQCK
      GO TO 729
  715 DVTA=0.2*(FGZ-FGKA)
      HQT=HQCK
  727 IF(0.0001-ABS(FGZ-FGKA))720,720,721
  720 HQT=HQT+DVTA
      CALL YPPKC(HQT,C1AO,HTO,B2O,AML1,ADST,
     *PIK,ALFA2,LA2,C2A,C2U,AT2,
     *II,RK,HA,RKI,HAI,STAGE,
     *AKP,XFK,XFA,UPR2,PI,PIT,ALFA4,LA4,C4A,HT,AKPX)
      FGZ=FGZA(ALFA2,ALFA4,LA2)
      IF(FGZ-FGKA)713,716,716
  713 B=1.
      A=0.5
      DVTA=-ABS(DVTA)*A
      GO TO 727
  716 IF(B.GT.0.5) A=0.5
      DVTA=ABS(DVTA)*A
      GO TO 727
  721 HQC=HQT
      GO TO 729
  712 CALL YPPKC(HQT,C1AO,HTO,B2O,AML1,ADST,
     *PIK,ALFA2,LA2,C2A,C2U,AT2,
     *II,RK,HA,RKI,HAI,STAGE,
     *AKP,XFK,XFA,UPR2,PI,PIT,ALFA4,LA4,C4A,HT,AKPX)
      FGZ=FGZA(ALFA2,ALFA4,LA2)
      IF(FGZ-FGKA)730,719,719
  730 ALIN=88
      HQC=88
      GO TO 729
  719 DVTA=0.2*(FGZ-FGKA)
  726 IF(0.0001-ABS(FGZ-FGKA))722,723,723
  722 HQT=HQT+DVTA
      CALL YPPKC(HQT,C1AO,HTO,B2O,AML1,ADST,
     *PIK,ALFA2,LA2,C2A,C2U,AT2,
     *II,RK,HA,RKI,HAI,STAGE,
     *AKP,XFK,XFA,UPR2,PI,PIT,ALFA4,LA4,C4A,HT,AKPX)
      FGZ=FGZA(ALFA2,ALFA4,LA2)
      IF(FGZ-FGKA)724,724,725
  724 B=1.
      A=0.5
      DVTA=-ABS(DVTA)*A
      GO TO 726
  725 IF(B.GT.0.5) A=0.5
      DVTA=ABS(DVTA)*A
      GO TO 726
  723 HQC=HQT
  729 CONTINUE
      IF(KRAN(22).EQ.0) GO TO 66
      WRITE(7,1357) HQC,ALIN
      WRITE(7,572)
   66 RETURN
      END
      ! 转子出口参数计算
      SUBROUTINE PPZPKC(HT,AKPX,QLA1,C1A,C2AOA,
     *PIK,ALFA2,LA2,C2A,C2U,AT2,
     *PI,PIT,ALFA4,LA4,C4A,II)       
      REAL LA2,LA4,LA2A
      COMMON/A2/R1/A6/R2,DK1,DK2,DK4,DDK1,DDK2,DDK4
      COMMON/A7/UK1/A8/ALFA1
      COMMON/A9/RRR,AT1,UK2,UK4/A10/AKH
      COMMON/A11/C2AO,C4AO,ALFA4O/A13/AK
      HZ=HT*UK1**2*AKH/9.81
      HF=HT*AKH
      C1U=C1A*COS(ALFA1)/SIN(ALFA1)  !进口周向相对速度系数
  305 CONTINUE
      DF2=DK2/DK1
      C2U=(HF+R1*C1U)/(R2*DF2**2)    !出口周向相对速度系数
      TAW=1-(C1U+C2U*DF2)/(2.0*R1)  !反动度
      IF(TAW-0.3)705,705,706
  705 TAW=0.3
      GO TO 708
  706 IF(0.8-TAW)707,707,708
  707 TAW=0.8
  708 TAW=TAW
      AKPDK=1-(1-AKPX)*(0.8*TAW+0.1)
      HAD=HZ*AKPX   !级的等熵功
      HADK=HZ*AKPDK  !转子的等熵功
      FPI=AK*AT1*RRR/(AK-1)
      FPP=AK/(AK-1)
      PI=(1+HAD/FPI)**FPP   !级压比
      PIK=(1+HADK/FPI)**FPP  !转子压比
      PIT=(1+HZ/FPI)**FPP    !级理论压比
      AT2=AT1+HZ/(FPP*RRR)   !出口总温
      SMUL=1.
      IF(HZ.LT.0.) SMUL=-1.
      IF(ABS(PI-1.).LT.1.E-5) PI=1.+SMUL*1.E-5
      IF(ABS(PI-1.).LT.1.E-5) PIK=1+SMUL*(HADK/HAD)*1.E-5
      IF(ABS(PI-1.).LT.1.E-5) PIT=1.+SMUL*(HZ/HAD)*1.E-5
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
      ! 理论压头系数的计算
      SUBROUTINE PHOY(HQ,UPR2,HTO,C1AO,B2O,H)
      COMMON/A1/D1/A2/R1/A5/HRK
      COMMON/A21/KRAN(40)
      BK(XX)=1.385-0.385*ALOG10(XX*10)
      BBK(XX,YY)=BK(XX)*(10*YY-5)+6-10*YY
      BBBK(XX,YY,ZZ)=BBK(XX,YY)*(3.5-5*ZZ)+BK(XX)*(5*ZZ-2.5)
      AKD(TT,PP)=0.6*TT/PP+0.4
      AKC(TT,PP)=0.3*TT/PP+0.7
  557 FORMAT('P H O Y 崁梹')
  558 FORMAT('P H O Y 妿崡厤')
 1357 FORMAT(6F20.6)
      HTC=HTO/C1AO
  101 IF(HQ-1)103,104,104
  103 IF(D1-0.5)105,105,106
  105 IF(HTC-0.5)107,107,108
  107 AKN=1.0
      GO TO 111
  108 IF(HTC-0.6)109,110,110
  109 AKN=BBK(HRK,HTC)
      GO TO 111
  110 AKN=BK(HRK)
  111 IF(UPR2-1)112,112,113
  112 AKE=AKD(HTC,0.5)
      GO TO 114
  113 AKE=AKC(HTC,0.6)
  114 IF(1-AKE)115,115,116
  115 AKE=1
  116 CONTINUE
      GO TO 140
  106 IF(D1-0.7)117,118,118
  117 IF(HTC-0.5)119,119,120
  119 AKN=(3.5-5*D1)+BK(HRK)*(5*D1-2.5)
      GO TO 123
  120 IF(HTC-0.6)121,122,122
  121 AKN=BBBK(HRK,HTC,D1)
      GO TO 123
  122 AKN=BK(HRK)
  123 IF(UPR2-1)124,124,125
  124 AKE=AKD(HTC,0.5)
      GO TO 126
  125 AKE=AKC(HTC,0.6)
  126 IF(1-AKE)127,127,128
  127 AKE=1
  128 CONTINUE
      AKE=AKE*(3.5-5*D1)+5*D1-2.5
      GO TO 140
  118 AKN=BK(HRK)
      AKE=1.0
  140 IF(UPR2-1)141,141,142
  141 DEL=0
      GO TO 143
  142 DEL=0.1
  143 H=(1-HQ)*AKN*AKE*(1+DEL-0.065*C1AO/(R1*(1-D1)*SIN(B2O)/COS(B2O)))
      GO TO 180
  104 IF(HQ-1.3)144,144,145
  144 H=0.9*(1-HQ)
      GO TO 180
  145 H=1.03-HQ
  180 CONTINUE
      IF(KRAN(10).EQ.0) GO TO 66
      WRITE(7,557)
      WRITE(7,1357) HQ,H
      WRITE(7,1357) AKN,AKE,DEL,C1AO,R1,D1,B2O
      WRITE(7,558)
   66 RETURN
      END
      FUNCTION PHMK(AKP,XFK)
      COMMON/A21/KRAN(40)
  559 FORMAT(14HP H M K 妿崡厤)
 1357 FORMAT(6F20.6)
      IF(0.5-XFK)187,188,188
  187 AKXF=1.5*XFK+0.25
      GO TO 1187
  188 AKXF=1
 1187 AKXF=AKP*AKXF
      PHMK=AKXF
      IF(KRAN(11).EQ.0) GO TO 66
      WRITE(7,1357) AKP,XFK,AKXF
      WRITE(7,559)
   66 RETURN
      END
      SUBROUTINE PHCY(HC,HQ,HQC,H)
      COMMON/A21/KRAN(40)
  560 FORMAT(14HP H C Y 妿崡厤)
 1357 FORMAT(6F20.6)
      HQT=HQ/HQC
      IF(HQT-0.4)149,150,150
  149 HK=1.362-1.8*HQT
      GO TO 151
  150 FF=5*HQT-2
      F=0.0595-0.007*FF
      F=F*FF
      F=F-0.3295
      F=F*FF
      HK=F+0.642
  151 H=HK+HQT*HC
      IF(KRAN(12).EQ.0) GO TO 66
      WRITE(7,1357) HC,HQC,H
      WRITE(7,560)
   66 RETURN
      END
      ! 效率偏离最佳点修正系数的确定
      SUBROUTINE PKPDOY(HQ,UPR2,AML1,HTO,VTA)
      COMMON/A10/AKH
      COMMON/A2/R1
      IF(HQ-1.0)43,44,44
   43   AML=AML1
      IF(UPR2-1.)45,45,46   !亚音与超音判定
   45 IF(AML1-0.5)47,48,48
   47 AML=0.5
   48 AP=4.55*AML-0.825
      BP=1.8*(AML-0.3)**3+0.08
      VTA=1-AP*(1.-HQ)**2-BP*(1.-HQ)
      GO TO 80
   46 IF(AML1-0.77)49,49,50
   49 BT=0.
      GO TO 39
   50 BT=1.47*AML1-1.1125
   39 VTA=1.-1.4*(1.-HQ)**2-BT*(1.-HQ)
      GO TO 80
   44 Y=((HQ-1.)*R1**2)/(HQ*HTO*AKH)
      IF(Y-0.35)51,51,52
   51 VTA=1.-0.5*Y**2
      GO TO 80
   52 IF(Y.GT..7) Y=.7
   53 VTA=1.-0.5*Y**2-2.2*(Y-0.35)**2
   80 CONTINUE
      RETURN
      END
      SUBROUTINE PKPDCY(HQT,VTAT)
      COMMON/A1/D1
      COMMON/A21/KRAN(40)
  564 FORMAT(18HP K P D C Y 妿崡厤)
 1357 FORMAT(6F20.6)
      IF(D1-0.55)57,58,58
   57 IF(HQT-0.6)59,60,60
   59 VTAT=0.916*HQT
      GO TO 67
   60 VTAT=1-1.1125*(1-HQT)
      GO TO 67
   58 IF(D1-0.73)61,62,62
   61 DD1=D1-0.55
      DDD=SQRT(DD1)
      IF(HQT-0.6)63,64,64
   63 VTAT=HQT*(0.916-0.28*DDD-0.5*DD1)
      GO TO 67
   64 AR=10.*DD1-8.5*DDD
       BR=4.*DDD-3.74*DD1+1.1125
      VTAT=1.-AR*(1.-HQT)**2-BR*(1.-HQT)
       GO TO 67
   62 IF(HQT-0.93)65,66,66
   65 VTAT=0.76*HQT
      GO TO 67
   66 VTAT=4.2*HQT-3.2
   67 CONTINUE
      IF(KRAN(16).EQ.0) GO TO 76
      WRITE(7,1357) D1
      WRITE(7,1357) HQT,VTAT
      WRITE(7,564)
   76 RETURN
      END
      ! 计算转子/静子级气动载荷的临界值
      SUBROUTINE FGCA(ANT,RE,FGK,FGKA)
      COMMON/A1/D1
      COMMON/A41/TTT
      COMMON/A21/KRAN(40)
  566 FORMAT(12H FGKA 崁墑厤)
 1357 FORMAT(6F20.6)
      AN=ANT          !相对换算转速
      RF=D1           !相对轮毂比
      IF(ANT.LT..74) ANT=.74
      IF(ANT.GT.1.1) ANT=1.1
      IF(D1.GT..7) D1=.7
      IF(D1.LT..4) D1=.4
      FN=2.089-1.269*ANT
      AF=(FN-1.3)/0.3
      BF=1.3-0.4*AF
      FGA=AF*D1+BF
      IF(FGA.GT.1.3) FGA=1.3
      IF(RE-163000)995,996,996
  995 FGR=1.3-0.184*RE/100000
      GO TO 999
  996 FGR=1
  999 FGK=FGA*FGR
      FGKA=FGK*1.2
      CONTINUE
      ANT=AN
      D1=RF
      IF(KRAN(18).EQ.0) GO TO 66
      WRITE(6,1357) FGK,FGKA
      WRITE(6,566)
   66 RETURN
      END
      ! 级的气动载荷临界值计算
      FUNCTION FGZA(ALFA2,ALFA4,LA2)
      REAL LA2
      COMMON/A14/HB3,HB4,CA,BT3,BT4
      COMMON/A21/KRAN(40)
  567 FORMAT(12HF G A 崁墑厤)
 1357 FORMAT(6F20.6)
      BF1=BT4+SIN(ALFA4)/HB4
      BF2=SIN(ALFA4)/BF1
      TF1=HB3*BT3+SIN(ALFA2)
      FF1=TF1*(BF2*TF1+HB3*SIN(ALFA2))
      FF2=TF1*BF2-HB3*SIN(ALFA2)
      FF3=(1.2-2*CA)/(1.0-0.35*LA2**2)
      FGZA=5.2*FF3*FF2/SQRT(FF1)
      FGZR=FGZA
      IF(KRAN(19).EQ.0) GO TO 66
      WRITE(7,1357) HB3,HB4,CA,BT3,BT4,ALFA2,ALFA4,LA2,
     *BF1,BF2,TF1,FF1,FF2,FF3
       WRITE(7,1357) FGZR
      WRITE(7,567)
   66 RETURN
      END
      ! 分离边界点相对流量系数的初值确定
      FUNCTION HQCP(C1AO,B2O,XXX,FGK)
      COMMON/IQ/IQP
      COMMON/A3/B,H1,H2,T1,T2,CK/A2/R1
      COMMON/A7/UK1/A8/ALFA1
      COMMON/A13/AK
      COMMON/A9/RRR,AT1,VK2,VK4
      BT2=B/T2
      BH2=B/H2
      BT1=B/T1
      BF=SIN(B2O)/(BT2+BH2*SIN(B2O))
      F=R1/C1AO-COS(ALFA1)/SIN(ALFA1)
      SBC=SQRT(1/(1+F**2))
      XX=XXX   !!!!!!!!!!!   k*R
      X1=SIN(ALFA1)
      X2=UK1**2/X1**2
      X3=(AK-1)/(2*XX)
      X4=X2*X3
      ATT=AT1-C1AO**2*X4
      XXX=SQRT(XX*ATT)
      AML1=C1AO*UK1/(SBC*XXX)
      AL1=SQRT((AK+1)*AML1**2/(2*(1+(AK-1)*AML1**2/2)))
  403 TF=H1*BT1+SBC
      IF(AL1.GE.0.9) AL1=0.9
      F1=1.0-0.35*AL1**2
      F2=6.24-10.4*CK  !CK转子最大相对厚度
      F3=TF*(BF*TF+H1*SBC)
      SBCC=(BF*TF-FGK*F1*SQRT(F3)/F2)/H1  !新的相对气流角正弦
	SB1=SBC
      SBC=SBCC
      F=SQRT(1/SBC**2-1)   !新的相对气流角余切值
      C1AC=R1/(F+COS(ALFA1)/SIN(ALFA1))  !流量系数
      IF(0.0001-ABS(SBC-SB1))401,402,402
  401 ATT=AT1-C1AC**2*X4
      XXX=SQRT(XX*ATT)
      AML1=C1AC*UK1/(SBC*XXX)
      AL1=SQRT((AK+1)*AML1**2/(2+(AK-1)*AML1**2))
      GO TO 403
  402 HQCK=C1AC/C1AO
      HQCP=HQCK
730   CONTINUE
      RETURN
      END
      SUBROUTINE YPPKC(HQ,C1AO,HTO,B2O,AML1,ADST,PIK,ALFA2,
     *LA2,C2A,C2U,AT2,II,RK,HA,RKI,HAI,STAGE,
     *AKP,XFK,XFA,UPR2,PI,PIT,ALFA4,LA4,C4A,HT,AKPX)
       REAL LA2,LA4
      COMMON/A2/R1/A6/R2,DK1,DK2,DK4,DDK1,DDK2,DDK4
      COMMON/A7/UK1/A8/ALFA1
      COMMON/A9/RRR,AT1,UK2,UK4/A10/AKH
      COMMON/A11/C2AO,C4AO,ALFA4O/A13/AK
      COMMON/A16/SIGMAK,SIGMAA,DSITAK,DSITAA
      COMMON/A21/KRAN(40)/R/INDEX5
      DIMENSION RK(20),HA(20),RKI(10),HAI(10),STAGE(24)
  570 FORMAT('Y P P K C 崁梹?')
  571 FORMAT(6F20.6)
      IF(KRAN(6).EQ.1) WRITE(7,570)
      IF(KRAN(6).EQ.1) WRITE(7,571) HQ,QLA1
      C1A=HQ*C1AO
      QLA1=QQL(C1A)
      IF(HQ.LE.1.) GO TO 1
      INDEX5=0
      CALL ROTOR(II,RK,HA,RKI,HAI,STAGE,HQ,QLA1,
     *C1U,C2A,C2U,C2,BETA2,ALFA2,LA2,AT2,ALFA4,LA4,C4A,
     *HZ,HZ1,HT,PIK,PIT,PI,AKPX,INDEX5)       ! 转子参数计算
      IF(INDEX5.EQ.1) RETURN
    1 CONTINUE
      CALL PHOY(HQ,UPR2,HTO,C1AO,B2O,H)
      H=H*PHMK(AKP,XFK)
      CALL PKPDOY(HQ,UPR2,AML1,HAI(10),VTA)   !主特性区效率计算
      CALL PKPDK(XFK,XFA,HQ,VTA)              !效率修正
      HT=H*R1**2+HQ*HTO
      AKPX=ADST*VTA
      IF(KRAN(6).EQ.1) WRITE(7,571) HT,AKPX,QLA1,C1A,UK4,HAI(9)
      C2AOA=HAI(9)
      CALL PPZPKC(HT,AKPX,QLA1,C1A,C2AOA,PIK,ALFA2,LA2,C2A,C2U,AT2,
     *PI,PIT,ALFA4,LA4,C4A,II)                   ! 静子出口参数计算
      RETURN
      END
      ! 根据叶型集合参数对效率修正系数的确定
      SUBROUTINE PKPDK(XFK,XFA,HQ,VTA)
      COMMON/A21/KRAN(40)
  563 FORMAT(16HP K P D K 妿崡厤)
 1357 FORMAT(6F20.6)
      IF(XFK-0.5)82,82,81
   81 AKVTA=20*XFK-9
      GO TO 97
   82 AKVTA=1
   97 IF(XFA-0.5)84,84,83
   83 AAVTA=20.*XFA-9.
      GO TO 85
   84 AAVTA=1
   85 AVTA=(AKVTA+AAVTA)*0.5
      IF(HQ-1)88,89,89
   88 VTA=1-(1-VTA)/AVTA
      GO TO 90
   89 VTA=1-(1-VTA)*AVTA
   90 CONTINUE
      IF(KRAN(15).EQ.0) GO TO 66
      WRITE(7,1357) HQ,XFK,XFA,VTA
      WRITE(7,563)
   66 RETURN
      END
      FUNCTION RLQMDA(GL2A)  !由折合流量求无因次速度数
      COMMON/A13/AK
      COMMON/A35/IJ
      IF(GL2A.GE.1) IJ=2
      IF(GL2A.GE.1.) GL2A=0.999999
      X=0.6*GL2A
      Y=1./(AK-1.)
    3 CONTINUE
      RLI=GL2A/(((AK+1.-(AK-1.)*X**2)/2.)**Y)
       Z=ABS(RLI-X)
      IF(Z-1E-6)2,2,1
    1 X=RLI
       GO TO 3
    2 RLQMDA=RLI
      RETURN
      END
      SUBROUTINE PCPBX
      DIMENSION D3(9),D31(30), D32(30)
      REAL ID1(33)
      DIMENSION DP(30,2),DC(30,2),DB(30,2),DBT(30,2),
     *R0(30,2),F0(30,2)
      EXTERNAL XHA,XHP,TXP
      EXTERNAL TXH,TKP
      COMMON/BB/D3,D31,D32,ID1
      COMMON/BB1/OB,PK,G3,PBX,T1,XTA3,AKP,IZ,  KF,KDB,UBX,PBHA
      COMMON/BB2/IZ1,IZ2,IZ3
      COMMON/Z8/DP,DC,DB
      COMMON/Z10/R0,F0,DBT
      COMMON/Z1/DP1,DC1,DB1,DPB,DCB,DBB/Z2/HZ1,HZC,HZB,HZK,
     *UK1,UKB/Z3/C1A,CAC,CAB/Z4/K,IV1
      COMMON/Z5/GP,PA,P,R/Z13/BC3
      COMMON/Z7/TAY,TAYC,TAYB
      COMMON/BX/ALF1,C1U,P1,CA,R1
    3 K1=0
      XA=XHA(T1)                 !进口总焓
      IF(K.LT.4)GOTO 11
      DP1=D31(1)
        DB1=D32(1)
        GOTO 28
   17 DP1=D3(1)                       !进口外径
        DB1=D3(3)                     !进口内径
      F1=FD(P,DP1,DB1)                !进口面积
        K1=3
        GOTO 24
   11 IF((D3(1).GT.0.).AND.(D3(3).GT.0.))GOTO 17
   13 C1=C1A/SIN(ALF1)                !进口绝对速度
        F1=F(C1A,C1,T1,P1,R,GP,XA)    !进口面积
       IF(K-2)14,34,6
   14 IF(K.EQ.0)GOTO 16
      DO 77  I=1,IZ2
      DP(I,1)=D31(I)
        DP1=D31(1)
   77 CONTINUE
      DO 19 I=1,IZ
      DP(I,2)=(D31(I)+D31(I+1))/2.
   19 CONTINUE
      GOTO 27
   34 DO 35 I=1,IZ2
      DC(I,1)=D31(I)
        DC1=D31(1)
   35 CONTINUE
      DO 36 I=1,IZ
      DC(I,2)=(D31(I)+D31(I+1))/2.
   36 CONTINUE
      GOTO 21
    6 DO 75 I=1,IZ2
      DB(I,1)=D31(I)
        DB1=D31(1)
   75 CONTINUE
      DO 7 I=1,IZ
      DB(I,2)=(D31(I)+D31(I+1))/2.
    7 CONTINUE
      GOTO 31
   16 DO 18 I=1,3
      IF(D3(I).GT.0.)IV=I
   18 CONTINUE
      IF(IV-2)25,23,30
   21 F12=-F1/2.
        DB1=DF(P,F12,DC1)
        K1=2
   22 DP1=DF(P,F1,DB1)
        GOTO (24,45),K1
   23 D1=D3(2)
        DP1=SQRT(4.*F1/P/(1.-D1**2))
         K1=2
        DB1=D1*DP1
   24 F12=F1/2.
      DC1=DF(P,F12,DB1)
      GOTO (45,46,32),K1
   25 DP1=D3(1)
   27 F12=-F1
        DB1=DF(P,F12,DP1)
   26 K1=1
        GOTO 24
   30 DB1=D3(3)
   31 K1=1
        GOTO 22
   28 DO 68 I=1,IZ2
        DP(I,1)=D31(I)
        DB(I,1)=D32(I)
      DBT(I,1)=DB(I,1)/DP(I,1)
      R0(I,1)=RD(DBT(I,1))
        DC(I,1)=2.*R0(I,1)*DP(I,1)
      F0(I,1)=FD(P,DP(I,1),DB(I,1))
   68 CONTINUE
      DO 69 I=1,IZ
        DB(I,2)=(D32(I)+D32(I+1))/2.
      DP(I,2)=(D31(I)+D31(I+1))/2.
      DC(I,2)=(DC(I,1)+DC(I+1,1))/2.
        R0(I,2)=(R0(I,1)+R0(I+1,1))/2.
      F0(I,2)=(F0(I,1)+F0(I+1,1))/2.
   69 CONTINUE
        F1=F0(1,1)
        DC1=DC(1,1)
   32 CONTINUE
      C1A=CAPOF(F1,T1,P1,R,GP,XA,ALF1)   !轴向速度
   45 D1=DB1/DP1
   46 R1=RD(D1)
   47 IF(K.EQ.4)GOTO 49
   48 DP(1,1)=DP1
        DC (1,1)=DC1
        DB(1,1)=DB1
        DBT(1,1)=D1
        R0(1,1)=R1
      F0(1,1)=F1
   49 IF(800.-Ob)52,52,54
   52 UK1=P*DP1*OB/60.
        GOTO 63
   54 UK1=OB
   63 C1U=C1A/TAN(ALF1)
        UC1=UK1*R1
      IF(KDB.GT.8)GOTO 58
      IF(IZ.EQ.1)GOTO 58
   62 HZ1=UK1**2*ID1(7)         !转子加功量
      C2U=C1U+HZ1/UC1           !出口绝对轴向速度
        TY=1.-(C1U+C2U)/2./UC1  !反动度
      TAY=TY
      GOTO 56
   55 C1U=UC1*(1.-TAY)-HZ1/2./UC1
      ALF1=ATAN 2(C1A,C1U)
      GOTO 3
   56 IF(TY-.5)57,58,58
   57 IF(ABS(TY-.5).LE..0001)GOTO 58
      TAY=.5
        GOTO 55
   58 CA=C1A/UK1
      RETURN
        END
      SUBROUTINE PCPBYX (PKP,PB,ALFB,AKGB,DBTB,TB,XA,XTA)
      DIMENSION D3(9),D31(30), D32(30)
      REAL ID1(33)
      DIMENSION DP(30,2),DC(30,2),DB(30,2),DBT(30,2),
     *R0(30,2),F0(30,2)
      COMMON/BB/D3,D31,D32,ID1
      COMMON/BB1/OB,PK,G3,PBX,T1,XTA3,AKP,IZ,  KF,KDB,UBX,PBHA
      COMMON/BB2/IZ1,IZ2,IZ3
      COMMON/Z8/Dp,DC,DB
      COMMON/Z10/R0,F0,DBT
      COMMON/Z5/GP,PA,P,R/Z13/BC3
      COMMON/Z1/DP1,DC1,DB1,DPB,DCB,DBB/Z2/HZ1,HZC,HZB,HZK,
     *UK1,UKB/Z3/C1A,CAC,CAB/Z4/K,IV1
      COMMON/Z7/TAY,TAYC,TAYB
      COMMON/Z17/FB
      COMMON/Z21/K10,K11,IDZ
   62 HAK=XHA(TXP(XHP(T1)+R*ALOG(PKP)))-XHA(T1)      !进出口焓升
   63 HZK=HAK/XTA
        XAB=XA+HZK                                   !出口总焓
        TB=TXH(XAB)
      IV1=0
      IW=0
   66 IF(K.LT.4)GOTO171
        DCB=DC(IZ2,1)
        DBTB=DBT(IZ2,1)
      DPB=D31(IZ2)
        DBB=D32(IZ2)
        FB=F0(IZ2,1)
   70 CAB=1.05*R*GP*TB/FB/PB
   71 CB=CAB/SIN(ALFB)
        FB1=F(CAB,CB,TB,PB,R,GP*AKGB,XAB)
      IF(K.EQ.4)GOTO 60
      IF(IW.EQ.3)GOTO 60
      GOTO 74
   60 IF(ABS(FB-FB1)-.0001)100,100,87
   87 IF(FB1-FB)72,100,187
  187 CAB=CAB*FB1/FB
        GOTO 71
   72 CAB=CAB*FB/FB1
        GOTO 71
  171 IF((D3(7).GT.0.).AND.(D3(9).GT.0.))GOTO 67
      GOTO 71
   64 IF((D3(7).GT.0.).AND.(D3(9).EQ.0.)) GO TO 65
      IF(ABS(D3(7)-D3(9)).LT..001)GOTO 103
      DBB=D3(9)
        IW=2
        GOTO 71
   65 DPB=D3(7)
        IW=1
        GOTO 71
   67 DPB=D3(7)
        DBB=D3(9)
        FB=FD(P,DPB,DBB)
      F1B=-FB/2.
        DCB=DF(P,F1B,DPB)
      DBTB=DBB/DPB
        IW=3
        GOTO 70
   74 FB=FB1
        IF(IW -1)79,86,77
   79 IF(K.EQ.0)GOTO 88
   75 IF(K-2)82,80,76
   76  DBB=D31(IZ2)
   77 DPB=DF(P,FB,DBB)
   78 F1B=FB/2.
        F1B=-F1B
      DCB=DF(P,F1B,DPB)
         GOTO 83
   80 DCB=D31(IZ2)
   81 F1B=FB/2.
        F1B=  -  F1B
      DBB=DF(P,F1B,DCB)
      DPB=DF(P,FB,DBB)
        GOTO 83
   82 DPB=D31(IZ2)
   86 F1B=-FB
        DBB=DF(P,F1B,DPB)
        GOTO 78
   88 DO 89 I=1,3
      IF(D3(I+3).GT.0.)IV1=I
   89 CONTINUE
      IF(IV1.EQ.0)GOTO 64
      IF(IZ1.EQ.0)GOTO (99,103,104),IV1
      GOTO (90,92,94),IV1
   90 DO 91 I=1,IZ3
        DO 91 J=1,2
   91 DP(I,J)=DP1
        GOTO 99
   92 DO 93 I=1,IZ3
        DO 93 J=1,2
   93 DC(I,J)=DC1
        GOTO 103
   94 DO 95 I=1,IZ3
        DO 95 J=1,2
   95 DB(I,J)=DB1
        GOTO 104
   99 DPB=DP1
        GOTO 86
  103 DCB=DC1
        GOTO 81
  104 DBB=DB1
        GOTO 77
   83 DBTB=DBB/DPB
        GOTO 100
  100 RETURN
        END
      SUBROUTINE ZC(IC)
      DIMENSION A(20,30),B(20,30),DP(30,2),DC(30,2),DB(30,2),
     *R0(30,2),F0(30,2),DBT(30,2),C1AS(30),C2AS(30)
      COMMON/MOD/IMODD(5),CAK(15),HZE(15),ALFK(15),CA2K(15),AKPK(15),
     *AKPDS(15)
      COMMON/Z1/DP1,DC1,DB1,DPB,DCB,DBB/Z2/HZ1,HZC,HZB,HZK,
     *UK1,UKB/Z3/C1A,CAC,CAB/CB/A,B/Z4/K,IV1
      COMMON/BB2/IZ1,IZ2,IZ3
      COMMON/Z7/TAY,TAYC,TAYB
      COMMON/Z8/DP,DC,DB
      COMMON/Z17/FB
      COMMON/Z21/K10,K11,IDZ
      COMMON/Z10/R0,F0,DBT
      COMMON/BX/ALF1,C1U,P1,CA,R1
      COMMON/C12AS/C1AS,C2AS
      UKC=(UK1+UKB)/2.
      M1=2
      HZZ=HZB
      HV1=HZ1
      CA1=C1A
      IF(IZ1.GT.0.AND.K.EQ.0)IC=IZ1
      IF((K.GT.0).OR.(IC.GT.0))GOTO 2
      IC=HZK/HZC/UKC**2+.5
    2 Z=IC
        Z0=Z-1.
        IC1=IC+1
        HZ0=HZK/Z
      IF(IC.EQ.1)GOTO 51
      IF(IC.EQ.2)GOTO 55
      IF(K11.GT.0)GOTO 48
   43  IF(IMODD(1).EQ.0)GOTO278
      HZZX=0.9*HZ0/(UKB**2)
      R1X=SQRT((1+(DB1/DP1)**2)/2.)
      RZX=SQRT((1+(DBB/DPB)**2)/2.)
      C1A1=C1A/UK1
      C1AZ=CAB/UKB
      HZ1X=HZZX*R1X*C1A1/RZX/C1AZ/1.1
      HZZ=HZZX*UKB**2
      HV1=HZ1X*UK1**2
278   CONTINUE
      IF(HZ0-(HV1+HZZ)/2)44,46,46
   44 HZZ=2.*HZ0-HV1
   46 B1=(HZK-Z*HV1-(HZZ-HV1)*Z*(2.*Z-1.)/6./Z0)/(Z*(Z-2.)/6./Z0)
      AHS=0.
      DO 47 I=1,IC
      AHS=AHS+A(20,I)
   47 CONTINUE
      DO 4 I=1,IC
      AI=I-1
      A(20,I)=A(20,I)*HZK/AHS
      IF(IMODD(1).GE.1.0.AND.IMODD(1).NE.10)A(20,I)=HZE(I)
    4 CONTINUE
   53 HZU=HZ0/UKC**2
        DP(IC1,1)=DPB
        DB(IC1,1)=DBB
        DC(IC1,1)=DCB
      DBT(IC1,1)=DBB/DPB
        R0(IC1,1)=RD(DBB/DPB)
        F0(IC1,1)=FB
      DO 7 I=1,IC1
        AI=I-1
      A(5,I)=2.*(CA1+CAB-2.*CAC)*(AI/Z)**2+
     *(4.*CAC-3.*CA1-CAB)
     **AI/Z+CA1
      A(5,I)=C1AS(I)
      IF(IMODD(2).EQ.1)A(5,I)=CAK(I)
      IF(K.EQ.0)GOTO 10
        IF(K-2)12,7 ,7
   10 IF(IV1-2)17,19,21
   12 A(1,I)=DP(I,1)
        A(2,I)=DP(I,2)
      B(1,I)=DP(I,2)
         B(2,I)=DP(I+1,1)
         GOTO 7
   17 IF(IV1.EQ.0)GOTO 26
      DO 18 J=1,2
   18 DP(I,J)=DP1
        GOTO 12
   19 DO 20 J=1,2
   20 DC(I,J)=DC1
        GOTO 7
   21 DO 23 J=1,2
   23 DB(I,J)=DB1
        GOTO 7
   26 DP(I,1)=.4*(DP1-DPB)*(AI/Z)**2-1.4*(DP1-DPB)*AI/Z+DP1
      GOTO 12
    7 CONTINUE
      GOTO 50
   54 A(20,1)=HZ1
        A(20,2)=HZB
        GOTO 53
   51 M1=1
      A(6,1)=CA1
       A(5,1)=CA1
       IF(HZ1.GT.0.)HZK=HZ1
        HZ1=HZK
        A(20,1)=HZ1
      GOTO 52
   48 Z9=IDZ
        Z6=Z-Z9
      HV1=HZ1*Z6/Z
        HZZ=HZB*Z6/Z
   52 UC1=UK1*DC1/DP1
      C2U=C1U+HV1/UC1
      TAY=1.-(C1U+C2U)/2./UC1
      IF(TAY.GE..5)GOTO 60
      IF(TAY.LT..4)GOTO 58
      TAY=.5
      C1U=UC1*(1.-TAY)-HV1/2./UC1
      ALF1=ATAN 2(C1A,C1U)
   60 GOTO (53,43),M1
   50 DO 5 I=1,IC
      IF(IC.EQ.1)GOTO 57
      AI=I-1
      A(19,I)=.34*(AI/Z0)**2-.41*AI/Z0+1.15
      B(19,I)=A(19,I)-.05
      IF(B(19,I).LT.1.)B(19,I)=1.
      A(8,I)=.93-.06*(AI/Z0)**2+.07*AI/Z0
      A(6,I)=(A(5,I)+A(5,I+1))/2.
      IF(IMODD(3).EQ.1)A(6,I)=CA2K(I)
   56 IF(IV1.EQ.0.AND.K.EQ.0)GOTO 77
      GOTO 5
   57 A(8,1)=.93
      IF(IC.GE.1)A(6,1)=(CA1+CAB)/2.
      IF(IMODD(3).EQ.1)A(6,1)=CA2K(1)
      A(19,1)=1.1
      B(19,1)=1.1
      GOTO 56
   77 DP(I,2)=(DP(I,1)+DP(I+1,1))/2.
              A(1,I)=DP(I,1)
        A(2,I)=DP(I,2)
      B(1,I)=DP(I,2)
        B(2,I)=DP(I+1,1)
    5 CONTINUE
      IF(IC.EQ.1)GOTO 58
      IF(IC.GT.2)GOTO 22
      A(9,2)=TAY+TAYB
      GOTO 58
   22 IF(IC-5)30,28,29
   30 N=2
        GOTO 31
   28 N=3
        GOTO 31
   29 IF(IC-10)32,32,33
   32 N=4
        GOTO 31
   33 N=5
        GOTO 31
   31 DO 35 I=1,N
        AN=N
         AI=I-1
      A(9,I)=TAY-AI*(TAY-TAYC)/(AN-1.)
   35 CONTINUE
   40 DO 37 J=N,IC
        AJ=J-N
      A(9,J)=TAYC+AJ* TAYB
      IF(A(9,J).GT..8)A(9,J)=.8
   37 CONTINUE
      IF(IMODD(4).NE.1)GOTO58
      DO581 JJ=2,IC
      A(9,JJ)=ALFK(JJ)
581   CONTINUE
   58 A(9,1)=TAY
      GOTO 38
   55 HZB=HZK-HZ1
      GOTO 54
  330 FORMAT(10X,12HTAU NIGE 0.5,2X,4HTAU=,F6.4/)
   38 RETURN
      END
      SUBROUTINE CAB5(KF,KC3,XF,B1,B2,X,CM,BT,
     *DB,B2T,DI,AI,AD,EP,J1,A,TT,B1K,B2K)
      PA=57.295828
      B11=(90.-B1)/100.
      SIB=SIN(B1/PA)
      B13=B1/100.
      AP=14.-7.*X
      AI=BT*(8.1*B11-.0659*B11**2+.308*B11**3-.53*B11**4)
      AD10=(11.5*B11**3/3-.25*B11**2+6.74*B11/3)*BT+(11.52*B11**3-
     *5.76*B11**4-6.04*B11**2+.28*B11)*BT**2
      AK=.02678*B11-.042447-1.2316*B11**2+2.2915*B11**3-
     *2.065*B11**4
      ALA=.1793+.39*B11+7.622*B11**2-20.05*B11**3+15.37*
     *B11**4
      AKY=.0612255*(4.4-BT)*(BT-1.)
      AN=AK+ALA*AKY
      AKD=142.*CM**3+9.34*CM**2+7.646*CM
      AMY=1.35*B13**3-3.125*B13**2+2.646*B13+.135
    7 AIC=18.217*CM-81.6861*CM**2-237.3737*CM**3+
     *2320.075*CM**4
      AIF=1.+1.6*(XF-.5)
      ADF=1.+.23*((2.*XF)**2-1.)/(.41-.002*B2)
    2 AM=.875*B13**2-.4375*B13**3-.6525*B13+.441
      IF(X.LE..775)GOTO 17
      AILA=39.359475*X-15.163398*X**2-21.396077
      ADLA=2.181818*X**2-.672727*X-.789091
      IF(AILA.GT.4.)AILA=4.
2929  FORMAT(1X,'2929',10F12.6)
      GOTO 18
   17 AILA=0.
      ADLA=0.
   18 IF(KC3.GT.0)GOTO 20
      IF(KF.EQ.0)GOTO 20
      IF(KF.EQ.3)FK=1.
    3 IF(KF.EQ.1)FK=1.05
      IF( KF.EQ.2)FK=1.1
    4 BTMY=AM/BT**AMY*ADF
      EP=(DB+1.0+ADLA-AILA-DI+FK*(AKD*
     *AD10-AIC*AI))/(1.-BTMY+AN*AIF)
      AI=FK*AIC*AI+EP*AN*AIF+DI+AILA-1.0
   21  AD=FK*AKD*AD10+BTMY*EP+ADLA
      IF(AI.LE.AP)GOTO 12
      AI=AP
      EP=(DB-AI+ADLA+FK*AKD*AD10)/(1.-BTMY)
      GOTO 21
   20 BTMY=AM/BT**AMY
      EP=(DB-AILA+ADLA+.7*(AKD*AD10-AIC*AI))/
     *(1.-BTMY+AN)
      AI=.7*AIC*AI+EP*AN+AILA
   22 AD=.7*AKD*AD10+BTMY*EP+ADLA
5050  FORMAT(1X,'5050',4F12.6)
      IF(AI.LE.AP)GOTO 12
      AI=AP
      EP=(DB-AI+ADLA+.7*AKD*AD10)/(1.-BTMY)
      GOTO 22
   12 B7=B2T
      B1K=B1+AI
      B2K=B7+AD
      EP=B2K-B1K
      TT=B1K+(1.5-2.*XF)*EP
      COT=COS(TT/PA)
      IF(KC3.GT.0)GOTO 14
   13 AK=.62-1.28*(COT/BT-.4)**2
      GOTO 15
   14 AK=.62-1.84*(COT/BT-.5)**2
      IF(J1.GT.0)GOTO 13
   15 E2=EP/2.
      SE2=SIN(E2/PA)
      STE=SIN((TT-E2)/PA)
      A=((SQRT(BT**2+4.*BT*STE*SE2+4.*SE2**2)
     *-BT)/2./SE2-AK*CM*BT)/SIB
      RETURN
      END
      SUBROUTINE KHXTA(G,H1,HD,HMI,   AG,KZ)
      DIMENSION G(4,30)
      DO 2 I=1,KZ
      H=H1-HD*(I-1)
      IF(H.LT.HMI)H=HMI
      G(3,I)=H
       G(2,I)=AG
    2 CONTINUE
      RETURN
      END
      SUBROUTINE BXHA(D,R,XF,BB,BT,H,A1,S,Q,ZA)
      P=3.14159
       PA=180./P
      D1=SQRT(2.*R**2-1.)
      Q=D*R
      Y=D*(1.-D1)/2.
      Z=P*BT*Q*H/Y
      I=Z+.5
      ZA=I
       HB=ZA/Z*H
       B=Y/HB
      AD=.25*(1.+(.23/(.41-.002*A1)*((2.*XF)**2-1.)))/BT**.965
      XPC=(90.-A1)/(.94+.025*BT-AD)
      AI=(.025*BT-.06)*XPC
      AD=AD*XPC
       ALK=A1-AD
      TET=ALK+(2.*XF-.5)*XPC
      A=1.+(BB-1.)*(R-D1)/(1.-D1)
      Q=B*(.9*BB/A*SIN(TET/PA)+S)
      BB=B
      RETURN
      END
      SUBROUTINE POFD(F,D1,D2,D3)
      COMMON/Z4/K,IV1
      P=3.14159
      IF(K.EQ.0.AND.IV1.EQ.0)GOTO 8
      IF((K.EQ.1).OR .(IV1.EQ.1))GOTO 8
      IF((K.EQ.2).OR.(IV1.EQ.2))GOTO 6
      IF((K.EQ.3).OR.(IV1.EQ.3))GOTO 9
    6 S=-F/2.
        D3=DF(P,S,D2)
        GOTO 7
    9 M=2
        GOTO 4
    8 M=1
        S=-F
        D3=DF(P,S,D1)
    4 S=F/2.
        D2=DF(P,S,D3)
        GOTO (15,7),M
    7 D1=DF(P,F,D3)
   15 RETURN
        END
      SUBROUTINE BAX(B2T,C21,B1,BAX2,BAX1)
      IF(C21.LT..85)C21=.85
      BAX2=ATAN(2./((C21-1.)/TAN(B1)/(1.+C21)+(3.*C21+1.)/
     *TAN(B2T)/(1.+C21)))
      BAX1=ATAN(2./((3.+C21)/TAN(B1)/(1.+C21)+(1.-C21)/
     *TAN(B2T)/(1.+C21)))
        RETURN
        END
      SUBROUTINE HIYX(IZ,KY)
      DIMENSION DP(30,2),DC(30,2),DB(30,2),DBT(30,2),
     *R0(30,2),F0(30,2),BT(30,2),BE1(30),BE2(30),
     *ALX(30),YD(30,2),IZKA(30,2)
      DIMENSION A(20,30),B(20,30)
      DIMENSION ZL(30,2),ZL0(30,2),BKA(30,2),YK(30)
      COMMON/CB/A,B
      COMMON/Z8/DP,DC,DB
      COMMON/Z10/R0,F0,DBT
      COMMON/Z11/BE1,BE2
      COMMON/Z12/BT
      COMMON/Z15/YD,ALX/Z16/IZKA,BKA
      COMMON/Z5/GP,PA,P,R/Z13/BC3
      COMMON/YKC/YK
      H1=YD(1,1)
      DH=H1-YD(IZ,1)
      DO 12 I=1,IZ
        Z0=IZ-1
        ZI=I-1
   11 IF(IZ.EQ.1)GOTO 9
      Z=ZI/Z0
      YKOH=H1-Z*DH
   13 IF(KY-2)1,3,15
   15 IF(KY-4)7,3,16
   16 IF(KY-6)7,5,7
    1 IF(IZ.EQ.1)AXM=5.
      IF(IZ.GT.1)
     *AXM=3.3*Z-3.*Z**2+5.3         !当量扩压器扩张角
    2 YAX=AXM**2*(BT(I,1)+COS(BE1(I)/PA)/2.)**2/2093.0625/((1.
     *+DBT(I,1))/R0(I,1))/BT(I,1)/(SQRT(F0(I,2)/F0(I,1)*
     *SIN(BE2(I)/PA))-SQRT(SIN(BE1(I)/PA)))**2    !展弦比
      IF(YAX.LE.YKOH)Y=YAX
        IF(YKOH.LE.YAX)Y=YKOH
    4 YD(I,1)=Y*YK(I)
    8 YD(I,2)=1*Y*YK(I)!原值1.15
      IF(I.EQ.IZ.AND.(B(14,IZ)-B(13,IZ)).GT.45.)
     *YD(I,2)=.6*Y
      GOTO 10
    3 IF(IZ.EQ.1)AXM=5.
      IF(IZ.GT.1)
     * AXM=3.1*Z-2.6*Z**2+5.3
      GOTO 2
    5 IF(IZ.EQ.1)AXM=5.
      IF(IZ.GT.1)
     *AXM=2.5*Z-2.2*Z**2+5.5
      GOTO 2
    7 IF(IZ.EQ.1)AXM=5.5
      IF(IZ.GT.1)
     *AXM=5.7+1.9*Z-1.8*Z**2
      GOTO 2
    9 YKOH=H1
      GOTO 13
   10 DO 12 J=1,2
      DC(I,J)=(DP(I,J)-DB(I,J))/2.    !叶高
      ZL(I,J)=P*BT(I,J)*DP(I,J)*R0(I,J)*YD(I,J)/DC(I,J)  !叶片数
      IZKA(I,J)=ZL(I,J)+.5
        ZL0(I,J)=IZKA(I,J)
      YD(I,J)=ZL0(I,J)/ZL(I,J)*YD(I,J)   !根据叶片数调整展弦比
      BKA(I,J)=DC(I,J)/YD(I,J)     !弦长
      A(9,I)=YD(I,1)
        B(9,I)=YD(I,2)
      IF(BKA(I,J).GT..015)GOTO 22
      BKA(I,J)=.015
        YD(I,J)=DC(I,J)/.015
   22 AXM=45.75*SQRT((1.+DBT(I,1))/R0(I,1)*BT(I,1)
     **YD(I,1))*(SQRT(F0(I,2)/F0(I,1)*
     *SIN(BE2(I)/PA))-SQRT(SIN(BE1(I)/PA)))
     */(BT(I,1)+COS(BE1(I)/PA)/2.)
      ALX(I)=AXM
   12 CONTINUE
      RETURN
        END
      SUBROUTINE PROH(J,A,DP,DB,CAI,CU,UK,CIG,PG,PW,R0)
      DIMENSION KP(30,3)
      DIMENSION A(20,30),DP(30,2),DB(30,2),CAI(30,2),
     *CU(30,2),UK(30),BE1(30),BE2(30),BT(30,2),PG(30,2),
     *BKA(30,2),IZKA(30,2),CIG(9),R0(30,2),PW(30,2)
      COMMON/Z11/BE1,BE2/Z12/BT/Z16/IZKA,BKA
      COMMON/PR1/KP
      CIG(1)=7800.
        CIG(2)=206E9
        CIG(3)=39E7
        CIG(4)=54E6
        CIG(5)=49E7
      DO 5 I=1,J
      D0=(DB(I,1)+DB(I,2))/(DP(I,1)+DP(I,2))
      D=DP(I+1,1)/DP(I,1)
        D=(D+3.)/4.
      D1=(1.-D0)*D*D/12
        C=D0*3+1.
      B1=(C+4.)*D1
        D1=D1*C
        B=CIG(1)*UK(I)**2
      X1=(CIG(3)/B-D1)/B1
        X2=(14.286*CIG(4)/B*
     *D0*D/BT(I,1)/R0(I,1)-D1)/B1
      D=(DP(I,1)+DP(I,2)-DB(I,1)-DB(I,2))**2/96
      D1=R0(I,1)
        B1=D1-D0
        C=D1*DP(I,1)/2
      D0=1.-D0
        D1=1.-D1
        C=6.2832*C/IZKA(I,1)
      B=(3.1*D0-3.15*B1)/D1*D*C
      D=(3.5*D0-3.75*B1)/D1*D*C
        C=PW(I,1)*CAI(I,1)
      D=D*((CAI(I,1)-CAI(I,2))*C-PG(I,2)+PG(I,1))
      B=B*C*(CU(I,1)-CU(I,2))
        C=0.0109083*(BE1(I)+BE2(I))
      C=B*SIN(C)+D*COS(C)
        D=(BKA(I,1)*0.1)**3
      C=C/D
    8 IF(X1-.2)1,2,2
    1 KP(I,1)=1
   2  IF(X2-.2)3,4,4
    3 KP(I,2)=1
    4 IF(CIG(5)-C)6,7,7
    6 KP(I,3)=1
    7 CONTINUE
    5 CONTINUE
      RETURN
        END
      SUBROUTINE ATAKC(D10,UK,BT,XF,CM,CA,CA2,ALA,A23,
     *AL1,BE1,BE2,B2T,BE11,BE21,EPS,AGA1,AI,DEL,TET,
     *Z,Z1)
      DATA PA/57.295828/
      IF(A23.GT.2.01)B2T=BE2
      DBE=BE2-BE1
        DBT=B2T-BE1
      R= SQRT((1.+D10**2)/2.)
      IF(ALA-.8)4,4,6
    4 X=(.4/XF)**2
      AKD=1.1
      AIL=AILA(R,UK,CA,AL1,BE1,ALA)
    2 DBE21=DBEO(B2T,BT)     !额定状态气流转角
    3 AI=BT*((BE1/30.)**2+1.5)+10.*D10/SQRT(BT)*(1.-DBT/
     *DBE21)-1.5-A23+AIL        !攻角
      IF(ALA.GT..8)AI=5.*(AI1-AI)*ALA+5.*AI-4.*AI1
      GOTO 7
    6 X=1.
        XF=.5
        AKD=1.
      GA2=ATAN  (1.2*CM)*PA
      AI=GA2+1.25
      IF(ALA.GT.1.)GOTO 7
      AI1=AI
      AIL=0.
      GOTO 2
    7 DC=(DBT+128.*X*CM)/(DBT+7.68*X)
      XM=.92*XF**2-.002*B2T+.18
      C21=CA2/CA
      X1=AKD*XM/SQRT(BT)*DCA(D10,C21)*DB2(B2T)*DC
      EPS=(DBT-AI)/(1.-X1)
      DEL=X1*EPS
      BE11=BE1+AI
      BE21=BE11+EPS
      TET=BE11+(1.5-2.*XF)*EPS
      COT=COS(TET/PA)
      IF(ALA-.85)11,11,10
   11 IF(BT.GT.1.45)GOTO 8
      IZ=Z*Z1+.5
      ZH=IZ
      Z1=ZH/Z
      BTH=BT*Z1
      Z=ZH
      BT=BTH
    8 AK=.62-1.28*(COT/BT-.4)**2
      GOTO 12
   10 AK=.62-1.84*(COT/BT-.5)**2
   12 EPS2=EPS/2.
      SEPS2=SIN(EPS2/PA)
      STE=SIN((TET-EPS2)/PA)
      AGA1=((SQRT(BT**2+4.*BT*STE*SEPS2+4.*SEPS2**2)-
     *BT)/2./SEPS2-AK*CM*BT)/SIN(BE1/PA)      !CM为Cmax
   20 FORMAT(8F8.4)
      RETURN
        END
      REAL FUNCTION DBEO(BE,BT)
      DBEO=(.001997*BE**2+.170433*BE+2.884)*(.3305+
     *.5*SQRT(3.2*BT-1.373))
      RETURN
        END
      REAL FUNCTION DCA(D1,C21)
      IF(D1-.4)2,2,4
    2 CA=1.
        GOTO 8
    4 IF(D1-.75)5,5,6
    5 CA=1.3*(D1-.4)/(C21-.545)-2.86*D1+2.145
      GOTO 8
    6 CA=.455/(C21-.545)
    8 IF(CA.GT.2.)CA=2.
      DCA=CA
       RETURN
      END
      REAL FUNCTION DB2(B2)
      IF(60.-B2)2,3,3
    2 DB2=1.
        GOTO 4
    3 DB2=1.+.0015*(60.-B2)+.0003*(60.-B2)**2
    4 RETURN
        END
      REAL FUNCTION AILA(R,UK,CA,ALF,BET,ALA)
      DATA PA/57.295828/
      IF(ALA-.4)2,2,4
    2 C=1.09
        GOTO 8
    4 IF(ALA-.81)5,5,6
    5 C=1.09-.375*(ALA-.4)**2-.075*(ALA-.4)
      GOTO 8
    6 C=2.53*ALA-1.93*ALA**2+.211
    8 AILA=ATAN  (1./(R*UK*C/CA-1./TAN(ALF   )))*PA-BET
   20 FORMAT(8F8.4)
      RETURN
        END
      SUBROUTINE PKIHA(I,U,P,T,C1,IZ,PBHA,CMKA)
      DIMENSION PG(30,2),PW(30,2),CIG(9)
      DIMENSION A(20,30),B(20,30),GKA(4,30),GB(10)
      DIMENSION DP(30,2),DC(30,2),DB(30,2),DBT(30,2),
     *R0(30,2),F0(30,2),CAI(30,2),BT(30,2),BE1(30),BE2(30),
     *AL2(30),AL4(30),DL(30,2),L(30,2),HZI(30),HAD(30),
     *PC(30,2),XT(30,2),AI(30,2),Q(30,2)
      REAL L,L1,L2
      DIMENSION CU(30,2),TE(30),AL1(30),UK(30),C1AS(30),C2AS(30)
      COMMON/C12AS/C1AS,C2AS
      COMMON/Z8/DP,DC,DB
      COMMON/Z10/R0,F0,DBT
      COMMON/Z11/BE1,BE2
      COMMON/Z12/BT
      COMMON/Z19/GKA,GB
      COMMON/PK1/ALFB
      COMMON/PK/CAI,AL2,AL4,AL1,DL,L,HZI,HAD,CU,TE,UK,PC,XT,Q,AI
      COMMON/Z1/DP1,DC1,DB1,DPB,DCB,DBB
      COMMON/CB/A,B/Z4/K,IV1
      COMMON/Z5/GP,PA,Q1,R/Z13/BC3
      COMMON/MOD/IMODD(5),CAK(15),HZE(15),ALFK(15),CA2K(15),AKPK(15),
     *AKPDC(15)
      COMMON/PP/PG,PW,CIG
      REAL*8 P,PC1,P1T,WORK
      D1=DBT(I,1)           !工作轮前
      R1=R0(I,1)
      ALF1=Al1(I)/PA
      IF(A(20,I).LE.1.)A(20,I)=A(20,I)*UK(I)**2
      HZ7=A(20,I)
      CA1=A(5,I)
      C2A=A(6,I)
      C4A=A(5,I+1)
      C =CA1/SIN(ALF1)
      B1=ATAN(CA1/(U-C1))
       BE1(I)=B1*PA
      AI1=XHA(T)
      AI(I,1)=AI1
      S1=XHP(T)
      AIT=AI1+.5*(U**2-2.*U*C1)
       WORK=EXP((S1-XHP(TXH(AI1-C**2/2.)))/R)
      PC1=P/WORK
      WORK=EXP((XHP(TXH(AIT))-S1)/R)
      P1T=P*WORK
      TKPT=TKP(TXH(AIT))
      AKP1=SQRT(2.*(AIT-XHA(TKPT  )))
      L1=CA1/SIN(B1)/AKP1
       L(I,1)=L1
      C2A=C2AS(I)
      A(6,I)=C2A
      CAI(I,2)=C2A
      W1=CA1/SIN(B1)
      CMK=.104-.06*L1
      CMK=CMK*CMKA
      IF(CMK.GT..08)CMK=.08
      AI2=HZ7+AI1             !工作轮后
       T2=TXH(AI2)            !出口总温
       S2=XHP(T2)             !出口熵函数
      C2U=(HZ7+C1*U)/U        !出口切向速度
      IF(IMODD(5).EQ.1)A(8,I)=AKPK(I)
      IF(IMODD(5).EQ.2)A(8,I)=(AKPDC(I)+1)/(2.0083+0.00007*(I-3))
      IF(IMODD(5).EQ.2)AKPK(I)=A(8,I)
      XTA=A(8,I)
  10  HA=HZ7*XTA
      AI2A=HA+AI1
      P2=P*EXP((XHP(TXH(AI2A))-S1)/R)     
      C2=SQRT(C2A**2+C2U**2)              !出口速度
        AIC2=AI2-C2**2/2.
      TC2=TXH(AIC2)                       !出口静温
      PC2=P2/EXP((S2-XHP(TC2))/R)         !出口静压
      IF(K.EQ.4)GOTO 25
 2    F2=F(C2A,C2,T2,P2,R,GP,AI2)         !出口面积
      CALL POFD(F2,DP(I,2),DC(I,2),DB(I,2))
      DP2=DP(I,2)
       DC2=DC(I,2)
       DB2=DB(I,2)
      D2=DB2/DP2
      R2=RD(D2)
       A(4,I)=D2
       B(3,I)=D2
      UK1=U/R1
      UC2=U*DC2/DC1
       UK2=UK1*DP2/DP1
      C2U=(HZ7+C1*U)/UC2
      B2=ATAN 2 (C2A,UC2-C2U)
        BE2(I)=B2*PA
      ABT=A(19,I)
        C21=C2A/CA1
        KC3=0
         A(10,I)=0.
         AK=.004
      AK=.0035
      IF(L1-.85)5,6,6
    6 KC3=1
       A(10,I)=1.
      V=I
      BC3I=BC3-.05*(V-1.)
      IF(BC3I.LT.1.)BC3I=1.
      GB(I)=BC3I
      BTK=(.8+3.*(HZ7/UK1**2-.1))/BC3I/R1*(1.+(BC3I-1.)*(R1-D1)
     */(1-D1))                   !?????????????????????????????????
    5 Z=IZ
      A(19,I)=ABT
      B(19,I)=ABT-.05
      IF(ABT.LT.1.)ABT=1
        ABT=1.1
    9 CONTINUE
      CALL BTIZI (BTK,C21,B1,B2,L1,CMK,DL(I,1),ZI,AK,ABT)
      IF(DL(I,1).GT..70)GOTO 24
      P2A=P*EXP((S2-S1)/R)                         !出口相对总压
      P2P=P2A*(1.-ZI*(P1T-PC1)/(P*EXP((XHP(TXH(
     *AIT+(UC2**2-U**2)/2.))-S1)/R)))              !出口绝对总压
      AITT=AIT+(UC2**2-U**2)/2.
      PC(I,1)=P2P/P
      S2A=S1+R*ALOG(PC(I,1))                       !出口等熵熵函数
      XTA=(XHA(TXP(S2A))-AI1)/HZ7                  !效率
      IF(IMODD(5).EQ.1)XTA=AKPK(I)
      IF(IMODD(5).EQ.2)XTA=AKPK(I)
      IF(ABS(A(8,I)-XTA).LT..0005)GOTO 12
      A(8,I)=XTA
      GOTO 10
   25 CONTINUE
   12 AL2(I)=ATAN 2 (C2A,C2U)*PA    
      ALF2=AL2(I)/PA
      KC3=0
      PG(I,1)=PC1
         PG(I,2)=PC2
      PW(I,1)=PC1/R/TXH(AI1-C**2/2)
      AKP2=SQRT(2.*(AI2-XHA(TKP(T2))))
      L2=C2/AKP2
      CMA=.14-L2*0.1
      CMA=CMA*CMKA
      IF(CMA.GT..08)CMA=.08
      BT(I,1)=BTK
         B(12,I)=CMA
         AI(I,2)=AI2
      CU(I,2)=C2U
         L(I,2)=L2
         F0(I,2)=F2
         R0(I,2)=R2
         DBT(I,2)=D2
      B(18,I)=P2P
      IS1=I+1                         !导向器后计算
      IF(I.NE.IZ)GOTO 13
      ALF4=ALFB
      C4U=C4A/TAN(ALF4)
      GOTO 14
   13 TAYB=A(9,I+1)
      C4U=UC2*(1.-TAYB)-A(20,I+1)/2./UC2
   14 IF(IMODD(5).EQ.0)XTAC=2.*XTA-1.
      IF(IMODD(5).EQ.2)XTAC=AKPDC(I)
   16 HAC=HZ7*XTAC
         AI4A=HAC+AI1
      P4=P*EXP((XHP(TXH(AI4A))-S1)/R)
      C4=SQRT(C4U**2+C4A**2)
      IF(I.EQ.IZ)GOTO 20
      F4=F(C4A,C4,T2,P4,R,GP,AI2)
      CALL POFD(F4,DP(IS1,1),DC(IS1,1),DB(IS1,1))
   20 DP4=DP(IS1,1)
         DC4=DC(IS1,1)
         DB4=DB(IS1,1)
      D4=DB4/DP4
         R4=RD(D4)
         A(4,I+1)=D4
         B(3,I+1)=D4
      UC4=U*DC4/DC1
         UK4=UK1*DP4/DP1
      IF(I.NE.IZ)GOTO 18
      C4U=C4A/TAN(ALF4)
      GOTO 19
18    IF(A(20,I+1).LE.1.)A(20,I+1)=A(20,I+1)*UK4**2
      C4U=UC4*(1.-TAYB)-A(20,I+1)/2./UC4
      ALF4=ATAN 2(C4A,C4U)
   19 C42=C4A/C2A
        AK=.005
      ABT=B(19,I)
      IF(ABT.LT.1.)ABT=1
        ABT=1.1
      CALL BTIZI(BTA,C42,ALF2,ALF4,L2,CMA,DL(I,2),ZIA,AK,ABT)
      P4=P2-ZIA*(P2-PC2)
      PC(I,2)=P4/P
      XT(I,2)=(XHA(TXP(S1+R*ALOG(PC(I,2))))-AI1)/HZ7
      IF(I.EQ.1)XT(1,2)=(XHA(TXP(S1+R*ALOG(PC(1,2)*PBHA)))-AI1)/HZ7
      IF(IMODD(5).EQ.2)XT(I,2)=AKPDC(I)
      IF(ABS(XT(I,2)-XTAC).LT..0005) GOTO 30
      XTAC=XT(I,2)
        GOTO 16
   30 CONTINUE
        HMM=(DP(I,1)-DB(I,1))/2
        DELX=L1**2*(0.025-HMM)/144.
        IF(HMM.GT.0.025)DELX=0.
        XTAC=XTAC-DELX
        XT(I,2)=XT(I,2)-DELX
      C2UT=(HZ7/GKA(3,I)+C1*U)/UC2
      PCU=UC2-C2UT
      B2T=ATAN 2 (C2A,PCU)
      IF(B2T.LT.0)B2T=B2T+3.14159
      CALL BAX(B2T,C21,B1,BX2,BX1)
         A(13,I)=BX1*PA
      IF(BX2*PA.LT.0.)BX2=(180.-BX2*PA)/PA
      A(14,I)=BX2*PA
        A(16,I)=(BX2-BX1)*PA
      IF(I.EQ.IZ)GOTO 21
      R0(I+1,1)=R4
        F0(I+1,1)=F4
        DBT(I+1,1)=D4
   21 CALL BAX(ALF4,C42,Alf2,AX4,AX2)
        B(13,I)=AX2*PA
      IF(AX4*PA.LT.0.)AX4=(180.+AX4*PA)/PA
      B(14,I)=AX4*PA
      Q(I,1)=C21
      Q(I,2)=C42
      UK(I+1)=UK4
      CU(I+1,1)=C4U
        AL1(I+1)=ALF4*PA
        A(18,I+1)=P4
        XT(I,1)=XTA
      AL4(I)=AL1(I+1)
      B(15,I)=B2T*PA
      B(1,I)=PC1
        B(2,I)=PC2
      BT(I,2)=BTA
        A(12,I)=CMK
      TE(I+1)=T2
      B(5,I)=zi
        B(6,I)=ziA
      GOTO 35
   24 GB(10)=2.
        GB(9)=I
   35 CONTINUE
      RETURN
        END
      FUNCTION QQL(C1A)
      REAL LA1
      COMMON/A15/AKR
      COMMON/A7/UK1/A8/ALFA1/A13/AK
      LA1=(C1A*UK1)/(AKR*SIN(ALFA1))
      FF=(AK+1)/2.
      FF1=1./(AK-1)
      FF2=(AK-1)/(AK+1)
      QQL=FF**FF1*(1-FF2*LA1**2.)**FF1*LA1
      QQ1=QQL
      RETURN
      END
      ! 级入口叶尖轴向相对速度的确定
      FUNCTION SPID(R1,BETA1,B1)
      SPID=R1/(CTG(BETA1)+CTG(B1))
      RETURN
      END
      SUBROUTINE AL1(C1A0,R1,B10S,AM10,AL10)
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
      FUNCTION CTG(F)
      CTG=1./TAN(F)
      RETURN
      END
      FUNCTION ARCCTG(F)
      IF(F.LT.0.000001) F=0.000001
      ARCCTG=ATAN(1./F)
      RETURN
      END
      FUNCTION TEMP(AT1)
      T=(AT1-200.)/1000.
      IF(AT1-600.)1,1,2
    1 TEMP=1.4015-0.156*T**2.-0.0022*T
      RETURN
    2 IF(AT1-800.)3,3,4
    3 TEMP=1.4202-0.111*T
      RETURN
    4 IF(AT1-1150.)5,5,6
    5 TEMP=1.4004-0.078*T
      RETURN
    6 TEMP=1.3686-0.0448*T
      RETURN
      END
      SUBROUTINE REYNOL(RK,C1AO,AKPDO,HTO,R1,AM10,RE,
     *C1AR,AKPDR,HTR)
      DIMENSION RK(20)
      AM=AM10
      X=ALOG(RE)
      IF(RK(16))3,4,3
    4 IF(AM10-.58)5,5,6
    5 AM10=.58
    6 CRE=.1977+.06365*X
      HRE=.7-.0556*X
      AKRE=(.195*AM10-.0589)*X+1.782-2.5*AM10
      HTR=HTO*CRE+HRE*R1**2.
      GOTO 10
    3 CRE=.0574*X+0.267
      IF(X-11.45)7,7,8
    7 HRE=.1*X-0.151
      GOTO9
    8 HRE=1.
    9 HTR=HTO*HRE
      AKRE=0.075*X+0.045
   10 C1AR=C1AO*CRE
      AKPDR=AKPDO*AKRE
      AM10=AM
      RETURN
      END
      FUNCTION Q(RWQ,QK)
      Q=((QK+1.)/2.)**(1./(QK-1.))*(1.-(QK-1.)*
     *RWQ**2/(QK+1.))**(1./(QK-1.))*RWQ
       RETURN
       END
      ! 最大入口轴向相对速度的确定
      SUBROUTINE BAHOB(PK,GDK,RIB,II)
      DIMENSION PK(20),GDK(6,30),RIB(4)
      COMMON/A7/UK1/A9/RRR,T1T,UK2,UK4/A13/AK/BLOK/P1T,A1,R,UKR
      COMMON/A22/C(4)
      COMMON/T1/T1
      AISP=0.
      RO=0.
       G=9.81
       CITP=0.995
    7 B1K=PK(5)
      B1T=PK(7)
      B=PK(8)
      HBX=PK(13)
      ALG=PK(14)
      R1=GDK(3,II)
      AGT=GDK(2,II)
      B1SP=GDK(6,II)      ! 进口相对气流角
      CAK=R1/(1./TAN(B1SP)+1./TAN(A1))
      CAM=1.03*CAK        ! 最大入口轴向相对速度初值
      X=2.*HBX            ! HBX叶片前缘厚度
      T=B/B1T
      AG=AGT*T
      AKP1=SQRT(2.*G*AK*R*T1T/(AK+1.))    ! 入口临界音速
      HSP=AG+HBX-T*SIN(B1SP)              ! 叶型前缘到叶栅喉部的距离
      DEH=-0.004542*X**4+0.04859*X**3-0.197*X**2+0.381*X
      D1=0.1
      EH=1.-HBX/(2.*T*SIN(B1SP))          ! 为入口轴向lambda数
      L=0
      IF(HSP.LT.0.) HSP=0.
      ICOUNT=0
   24 ICOUNT=ICOUNT+1
      AISP=0.5*ATAN(2.*ALG*HSP/(ALG**2-HSP**2))   ! 相对几何冲角初值
      CALL AL1(CAM,R1,B1,AM1S,RWQ1S)          ! 入口气流的马赫数及速度系数的确定
       IF(L.GT.3) GO TO 23
      W1=CAM*UK1/SIN(B1)                      ! 入口气流相对速度
      AGA1=AGT/SIN(B1)
      RWQ1=CAM*UK1/SIN(A1)/AKP1               ! RWQ1:进口绝对速度系数  RWQ1S：相对值
      AMU=AMQ(T1)
      RO=P1T*(1.-(AK-1.)*RWQ1**2/(AK+1.))**(1./(AK-1.))/(G*R*T1T) ! 入口气流密度
      RE=RO*W1*ALG*10.**(-3.)/AMU             ! 入口雷诺数
      DESP=0.04625*ALG/RE**0.2
      AIH=ASIN((0.9+0.0928*1.151**(10.*RWQ1S-9.)*(101.-100.*EH)**
     *(0.9153*RWQ1S-0.6985))*SIN(B1))-B1      ! 计算攻角的经验公式
      IF(AIH.LE.0.) AIH=0.
      AIBX=B1K-B1                             ! 设计点冲角，由己知几何数据和气动参数计算得到的攻角
      IF(AIBX.GT.AIH) AIH=AIBX
      AIC=ABS(AISP)-ABS(AIBX)                 ! 冲角变化量
      IF(AIBX.GE.0.) AIC=AIH+AISP
      ZWQ1S=0.5*(RWQ1S+1./RWQ1S)              ! 冲量函数
      DEO=(DESP+DEH)/AG
      ZWQPS=AK*RWQ1S*COS(AIC)/(AK+1.)+(ZWQ1S-AK*RWQ1S/(AK+1.))*
     *SIN(B1+AIC)/SIN(B1)                     ! 据冲角变化修正的冲量函数
      IF(ZWQPS.LT.1) ZWQPS=1.
      RWQPS=ZWQPS-SQRT(ZWQPS**2-1.)           ! 入口轴向速度系数修正值
      QWQ1S=Q(RWQ1S,AK)
      QWQPS=Q(RWQPS,AK)
      X3=10.*((1.-QWQPS)/QWQPS)
      CIP=QWQ1S*SIN(B1)/(QWQPS*SIN(B1+AIC))
      EWQPS=(1.-(AK-1.)*RWQPS**2/(AK+1.))**(1./(AK-1.))           ! 密度函数
      DZ=0.00321*X3**4-0.03775*X3**3+0.1743*X3**2-0.4138*X3+0.5   ! 修正轴向相对速度的经验系数
      IF(X3.GT.4.) DZ=0.039
      CIR=1.-DZ*(1.-QWQPS)*AK*EWQPS/(AK+1.)
      AGEOM=CITP*CIP*CIR*(1.-DEO)
      AZAP=QWQ1S/AGA1
      D2=AGEOM-AZAP           ! 最大轴向相对速度的修正值
      D3=ABS(D2)-ABS(D1)
      IF(ABS(D2).LT.ABS(D1)) GO TO 17
   19 D1=0.5*D2
      GO TO 18
   17 D1=D2
      IF(ABS(D3).LT.0.002) GO TO 19
   18 CONTINUE
      IF(ABS(D1).LT.0.001) GO TO 22   ! 精度满足要求，计算结束
      XCAM=CAM+D1                     ! 修正后的最大轴向相对速度
      IF(XCAM.GT.CAM*1.1) XCAM=CAM*1.1
      CAM=XCAM
      IF(ICOUNT.EQ.20) GOTO 22        ! 循环次数超过20次，计算结束
      GOTO 24
      ! 计算入口临界速度
   22 AM1A=AM1S*SIN(B1)   !绝对马赫数
      DCA=3.55*AM1A-2.65*AM1A**2-0.189
      IF(AM1A.LT.0.67) GO TO 23
      CAM=CAM*DCA   !轴向马赫数对最大流量系数修正
      L=9
      GO TO 24
   23 IF(RWQ1S.LT.0.3) RWQ1S=0.3
      A=3.*RWQ1S+1.78
      B0=1.+0.047/RWQ1S**A                    !计算缩放因子
      IF(B0.LT.1.02) B0=1.02
      B1O=ASIN(SIN(B1)/B0)-PK(18)/57.296      ! 最优入口相对气流角
      C1AOK=R1/(1./TAN(B1O)+1./TAN(A1))       !最优流量系数
      AA=3.82*RWQ1S+1.434   
      BKP=1.+0.017/RWQ1S**AA                  !缩放因子
      IF(BKP.LT.1.02) BKP=1.02
      AIOK=B1K -B1O                           ! 最佳攻角
      B1KP=ASIN(SIN(B1)/BKP)                  ! 经验确定临界入口相对气流角
      C1AKPK=R1/(1./TAN(B1KP)+1./TAN(A1))     ! 入口临界相对轴向速度
      RIB(1)=CAM
      RIB(2)=C1AKPK
      RIB(3)=AIOK
      RIB(4)=C1AOK
      C(1)=RIB(1)
      C(2)=RIB(2)
      C(3)=B1
      C(4)=B1KP
      RETURN
      END
      ! 计算理论能头
      SUBROUTINE HT(RK,GK,R1,R2,C1A0,C2C1,AL10,BETA2,ADST,
     *HT0,TAW,AT2,C2U,PI,PIK,HAD,HADK,HZ,C2A0,AL20,AM20,
     *PIT0,SIGMAK,SIGMAA,DSITAK,DSITAA,AKP2,GL10,HTS,IFAKT)
      DIMENSION RK(20),GK(4)
      COMMON/A7/UK1/A9/RRR,AT1,UK2,UK4/A13/AK
      COMMON/BLOK/AP1,A1,R,UKR
      HT0=R2**2*(RK(2)/RK(1))**2-C1A0*(R1*CTG(A1)+R2*(RK(2)/RK(1))
     ***2*C2C1*CTG(BETA2))                            !理论能头
      HT0=HT0*HTS
      C1U=C1A0*CTG(A1)                                ! 级进出口相对周向速度
      C2U=(HT0*GK(3)+R1*C1U)/(R2*(RK(2)/RK(1))**2)    ! 级进出口相对周向速度
      TAW=1.-(C1U+C2U*RK(2)/RK(1))/(2*R1)             ! 确定转子和静子效率的一个修正系数
      HZ0=HT0*UK1**2*GK(3)/9.81                       ! 理论压头校正系数
      AT2=AT1+HZ0*(AK-1)/(AK*R)                       ! 出口总温
      IF(TAW-0.8)1,2,2
    2 TAW=.8
      GOTO4
    1 IF(TAW-0.3)3,3,4
    3 TAW=.3
    4 AKPDM=1.-(1.-ADST)*(0.8*TAW+0.1)                ! 级效率最优时转子效率
      HAD=HZ0*ADST                                    ! 级等熵功
      HADK=HZ0*AKPDM                                  ! 转子等熵功（级加功量）
      PI=(HAD*(AK-1.)/(AK*R*AT1)+1.)**(AK/(AK-1.))    ! 级压比
      PIK=(HADK*(AK-1.)/(AK*R*AT1)+1.)**(AK/(AK-1.))  ! 转子压比
      T2A=AT2-(C2U*UK2)**2.*(AK-1.)/(2.*AK*9.81*R)    ! 出口轴向总温
      AKP1=SQRT(2.*9.81*AK*R*AT1/(AK+1.))             !进口临界速速
      XL10=C1A0*UK1/(SIN(A1)*AKP1)                    !进口切向无因次速度数
      GL10=XL10*((1-(AK-1)*XL10**2./(AK+1))*(AK+1)/2.)
     ***(1./(AK-1.))                                  !转子进口折合流量
      GL2A=GL10*SIN(A1)*(RK(1)/RK(2))**2.*((1-RK(3)**2.)/
     *(1-RK(4)**2.))*SQRT(AT2/AT1)*((AT2/T2A)**((AK+1)/
     *(2.*(AK-1))))/PIK                               !转子出口轴向折合流量
      IF(GL2A.GE.1) IFAKT=1   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF(GL2A.GE.1) GL2A=0.999999
      XL2A=RLQMdA(GL2A)                               !出口轴向无因次速度数
      C2A0=XL2A*SQRT(2.*9.81*AK*R*T2A/(AK+1))/UK2
      AKP2=SQRT(2.*9.81*AK*R*AT2/(AK+1))
      AL20=UK2*SQRT(C2A0**2.+C2U**2.)/AKP2            !出口无因次速度数
      AM20=AL20*SQRT(2./(AK+1))/SQRT(1-((AK-1)*AL20**2.)/(AK+1))  !马赫数
      PIT0=(1.+HZ0*(AK-1)/(AK*R*AT1))**(AK/(AK-1))  !??????级等熵压比
      IF(HZ0-10)7,8,8
    7 SIGMAK=1.+(HADK-HZ0)/(R*AT1)
      SIGMAA=1.+(HAD-HADK)/(R*AT1)
      GOTO9
    8 SIGMAK=PIK/PIT0   !转子总压恢复系数
      SIGMAA=PI/PIK     !静子总压恢复系数
    9 EPSIL1=(1.-((AK-1.)/(AK+1.))*AL10**2.)**(1./(AK-1.))
      EPSIL2=(1.-((AK-1.)/(AK+1.))*AL20**2.)**(1./(AK-1.))
      DSITAK=(1.-SIGMAK)*(AK+1)/((AL10**2.)*AK*EPSIL1)   !转子总压损失系数
      DSITAA=(1.-SIGMAA)*(AK+1)/((AL20**2.)*AK*EPSIL2)
      HZ=HZ0
      RETURN
      END
      SUBROUTINE POINT(RK,HA,GK,GB,GDK,
     *GDA,HAI,RKI,STAGE,II,INDEX1,OPTIMS,M)       ! 
      DIMENSION OPTIMS(10)
       DIMENSION RK(20),HA(20),GK(4),GB(10),GDK(6,30),
     *GDA(6,30),STAGE(24),RKI(10),HAI(10),RIB(4)
      COMMON/A16/SIGMAK,SIGMAA,DSITAK,DSITAA
      COMMON/A7/UK1/A9/RRR,AT1,UK2,UK4/A13/AK
      COMMON/A20/ARKI(30,10),AHAI(30,11)
      COMMON/A22/C(4)
      COMMON/BLOK/AP1,A1,R,UKR
      COMMON/T1/T1
   71 FORMAT(6F12.4)
  151 FORMAT('   崊?憰巹垖帒拡 帍拡寑嫓崨?弨悁寘拹巶')
      GB(4)=GB(4)/57.296
      IFAKT=0          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      AHP=RK(1)*(1.-RK(3))/2.         !叶高
      SK1=HA(13)/AHP*100              !径向间隙占叶高百分比
      SK1=SK1-1.
      IF(SK1.LE.0.) SK1=0.
      VTAC=0.444/(SK1+3.333)+0.867    !径向间隙对效率的影响
      HTS=1.4/(SK1+5.)+0.72           !径向间隙对加功量的影响
      AP1=GB(1)
      AT1=GB(2)
      UKR=GB(3)
      A1=GB(4)                        !进口绝对气流角
      R=GB(6)
      EK=GDK(1,II)
      EA=GDA(1,II)
      R1=GDK(3,II)
      R2=GDK(4,II)
      R4=GDA(4,II)
      B1=GB(4)
    2 DO 4 J=1,50
      IF(INDEX1-1)308,306,306
  306 IF(J-2)4,307,6
  308 IF(J-1)5,5,6
  307 C2C1=AHAI(II,1)/ARKI(II,1)
      C4C2=AHAI(II,11)/AHAI(II,1)
      DO 302 L=1,10
      RKI(L)=ARKI(II,L)   !?????????????
      HAI(L)=AHAI(II,L)
  302 CONTINUE
      GO TO 6
    5 B20=0.
      A40=0.
      C2C1=1.
      C4C2=1.
      GOTO7
    6 B20=RKI(2)
      A40=HAI(2)
    7 CALL DELTAK(RK,GK,J,EK,C2C1,B20,DELTK,BETA2)    ! 转子落后角计算
      RKI(2)=BETA2
      CALL DELTAK(HA,GK,J,EA,C4C2,A40,DELTA,ALFA4)    ! 静子落后角计算
      HAI(2)=ALFA4
      IF(RK(16))3,8,3
    8 IF(J-2)9,9 ,11
    9 CONTINUE
      WY=0.
      CALL AI0(RK,WY,BETA2,QIK,BETA1,2.)              ! 冲角
      QK=SPID(R1,BETA1,B1)
      IF(INDEX1-1)12,11,11         !????????
   11 DIK=RKI(8)
      CALL AI0(RK,DIK,BETA2,AI0K,BETA1,2.)
      RKI(4)=BETA1                !进口气流角
   12 C1A0K=SPID(R1,BETA1,B1)     !转子最优点进口流量系数
      GOTO40
    3  CONTINUE
      CALL BAHOB(RK,GDK,RIB,II)   ! 最大入口轴向相对速度的确定
      C1AM=RIB(1)
      C1AKR=RIB(2)
      AI0K=RIB(3)
      C1A0K=RIB(4)
      BETA1=RK(5)-AI0K
   40 IF(J-2)41,41,43
   41 CONTINUE
        US=0.
        YS=3.
        CALL AI0(HA,US,ALFA4,QIA,ALFA3,YS)
      QA=SPID(R2,BETA2,ALFA3)
      IF(INDEX1-1)44,43,43
   43 DIA=HAI(8)
        YS=3
      CALL AI0(HA,DIA,ALFA4,AI0A,ALFA3,YS)
      HAI(3)=AI0A
      HAI(4)=ALFA3
   44 C2A0A=SPID(R2,BETA2,ALFA3)   !静子最优点进口流量系数
      IF(J-1)13,13,14
   13 C1A0=C1A0K
      GO TO 16
   14 CONTINUE
      C1A0=RKI(1)
   16 CONTINUE
      CALL AL1(C1A0,R1,B10S,AM10,AL10)    !级最优点相对速度系数
      AKPDK=AKPD(RK,AL10)                 !转子效率
      HLIM=HP(AHP,AL10)                   !叶高对效率修正系数
      IF(J-1)19,19,20
   19 ADST=(AKPDK-HLIM)*GK(4)*VTAC        !最大效率校正系数
      GO TO 22
   20 ADST=RKI(10)
   22 CONTINUE
      CALL HT(RK,GK,R1,R2,C1A0,C2C1,AL10,BETA2,ADST,
     *HT0,TAW,AT2,C2U,PI,PIK,HAD,HADK,HZ,C2A0,AL20,AM20,
     *PIT0,SIGMAK,SIGMAA,DSITAK,DSITAA,AKP2,GL10,HTS,IFAKT)   ! 计算理论能头
      HAI(6)=AM20
      HAI(1)=C2A0
      HAI(10)=HT0
      C2C1=C2A0/C1A0
      IF(RK(16))45,46,45
   46 RKI(8)=ATAKA(RK,QIK,R1,QK,B1,AL10)    !级最优点速度系数对冲角修正
   45 HAI(8)=ATAKA(HA,QIA,R2,QA,BETA2,AL20)
      CALL OKPD(RK,C1A0,C1A0K,AL10,ANWK,OK,FK)  !转子效率修正系数
      CALL OKPD(HA,C2A0,C2A0A,AL20,ANWA,OA,FA)
      RKI(9)=C1A0K
      HAI(9)=C2A0A
      RKI(7)=AL10
      HAI(7)=AL20
      AKPDA=AKPD(HA,AL20)      !静子效率              VTAC径向间隙对效率修正
      RKI(10)=((AKPDK*OK*TAW+AKPDA*OA*(1.-TAW))-HLIM)*GK(4)*VTAC !级最优效率
      OC1A0=ONWK(RK,ANWK,ANWA,TAW,AKPDK,AKPDA,FK,FA) !级与转子最优点流量系数比
      C1A0=C1A0K*OC1A0    !级最优点流量系数
      IF(J.EQ.1)RKI(1)=C1A0
      IF(RK(16).EQ.0.) GO TO 26
      IF(C1A0.GT.C1AM) C1A0=C1AM
   26 CONTINUE
      RKI(1)=C1A0
      GL4A=GL10*SIN (A1)*(RK(1)/HA(2))**2.*(1-RK(3)**2.)*
     *SQRT(AT2/AT1)/(SIN(ALFA4)*(1-HA(4)**2.)*PI)    !级（静子）出口流量系数
      IF(GL4A.GE.1) IFAKT=1    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF(GL4A.GE.1) GL4A=0.999999
      AL4A=RLQMDA(GL4A)
      C4A0=AL4A*AKP2*SIN(ALFA4)/UK4   !级（静子）出口流量系数
      C4C2=C4A0/C2A0
      SA1=SIN(A1)
      SA4=SIN(ALFA4)
      ST21=SQRT(AT2/AT1)
  500 CONTINUE
      STAGE(19)=C4A0
      IF(HT0.LT.0.07) IFAKT=1  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF(IFAKT.EQ.1.AND.J.GT.2) GK(1)=GK(1)-0.1
      IF(IFAKT.EQ.1.AND.J.GT.2) OPTIMS(II)=II
      IFAKT=0   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      IF(J-2)30,30,28
   28 X=ABS(RKI(1)-STAGE(1))
      Y=ABS(RKI(10)-STAGE(4))
      Z=ABS(HAI(10)-STAGE(5))
      IF(X-0.001)29,29,30
   29 IF(Y-0.001)31,31,30
   31 IF(Z-0.001)32,32,30
   30 STAGE(1)=RKI(1)
      STAGE(4)=RKI(10)
    4 STAGE(5)=HAI(10)
      PRINT 151
      PRINT 71,(RKI(IS),IS=1,10)
      PRINT 71,(HAI(IS),IS=1,10)
      PRINT 71,X,Y,Z
      OPTIMS(M)=100*II
   32 B=SQRT(2.*AK/(9.81*R*(AK+1)))*1./(10.**3.*((AK+1)/2.)**(1./
     *(AK-1)))
      AMW=AMQ(T1)   !动力粘度
      RE=B*AP1*RK(8)*GL10*SIN(A1)/(AMW*SQRT(AT1)*SIN(BETA1)) !雷诺数
      IF(HA(17).LT.0.5) RE=1E6
      IF(RE-3.5*10.**5)33,33,34
   33 C1AO=RKI(1)
      AKPDO=RKI(10)
      HTO=HAI(10)
      CALL REYNOL(RK,C1AO,AKPDO,HTO,R1,AM10,RE,C1AR,AKPDR,HTR)   !雷诺数影响
      STAGE(1)=C1AR
      STAGE(4)=AKPDR
      STAGE(5)=HTR
      STAGE(23)=C1AR/C1AO
      GOTO35
   34 STAGE(1)=RKI(1)
      STAGE(4)=RKI(10)
      STAGE(5)=HAI(10)
      STAGE(23)=1.
   35 STAGE(2)=HAI(2)
      STAGE(3)=AT2
      STAGE(6)=AM10
      STAGE(7)=PI
      STAGE(8)=PIK
      STAGE(9)=HAI(1)
      STAGE(10)=C2U
      STAGE(11)=HAD
      STAGE(12)=HADK
      STAGE(13)=HZ
      STAGE(14)=C1A0K
      STAGE(15)=RKI(2)
      STAGE(16)=AI0K
      STAGE(18)=AK
      STAGE(20)=DELTK
      STAGE(21)=DELTA
      STAGE(17)=RE
      RKI(3)=STAGE(16)
      RKI(5)=STAGE(3)
      RKI(6)=STAGE(6)
      HAI(5)=STAGE(3)
      RETURN
      END
      SUBROUTINE BTIZI(BT,CC,B1,B2,L,CM,D,Z,AK,ABT)
      REAL L
      COMMON/MOD/IMODD(5),CAK(15),HZE(15),ALFK(15),CA2K(15),AKPK(15),
     *AKPDC(15)
      COMMON/ZIT/IZI,ALN,ALW,PKN,PKW
      REAL*8 AKP,Z1,AL
      S1=SIN(B1)
        S2=SIN(B2)
      IF(L.GE.1.)GOTO 3
      A=.634-CC*S1/S2
        BTT=(A+SQRT(A**2+.183*(COS(B1)*S2-CC*S1
     **COS(B2))/S2))/.183
        BD=BTT*ABT
    1 IF(BD.LT.1.)BD=1.
        IF(BD.GT.1.8)BD=1.8
      IF(AK.GT..0045)GOTO 17
      IF(L.GE..85)GOTO 16
   17 BT=BD
        GOTO 3
   16 BT=BD+6.66667*(L-.85)*(BT-BD)
3     D=1.-CC*S1/S2+S1*(1/TAN(B1)-CC/TAN(B2))/2./BT
      IF(D.LT..05)D=.05
      AKP=.0625*D**3-.025*D**2+.0225*D+AK
      IF(BT-1.)2,4,4
    2 X=.25
        GOTO 5
    4 X=.7
    5 Z1=AKP*2.*BT**X/S2**1.2
      IF(L-.5345)7,8,8
    7 AL=1.
        GOTO 10
    8 AL=1.+1.35*D**(.85/L)*(2.*L-1.)
   10 CONTINUE
      ACM=1.+L*(CM*10.)**5.68
      Z=AL*ACM*Z1
        IF(IZI.EQ.0)GOTO132
        PK=PKW+(PKN-PKW)*(ALW-L1)/(ALW-ALN)
        IF(L1.GT.ALW)PK=PKW
        IF(L1.LT.ALN)PK=PKN
        Z=Z*PK
132     CONTINUE
      RETURN
        END
      SUBROUTINE BRIDGE(II,STAGE,FX,GB)
      COMMON/A8/ALFA1
      COMMON/A9/RRR,AT1,UK2,UK4
      COMMON/A23/INDEX6,INDEX7,INDEX8,INDEX9,GUIDE
      COMMON/A24/IL,IM,IN,KM,KMI
      COMMON/A25/SL1(7),SL2(6),DL1(11),DL2(10)
      COMMON/A27R/SIGKKP,SGKP,SMKP
      COMMON/A25R/PDL1(11),PDL2(10)
      COMMON/A27/SIGMA,DSIGMA,DELTQ
      COMMON/A28/BET2,PPI
      COMMON/A29/AT2
      COMMON/A34/PM
      COMMON/C2/KSI,SL21(6)
      DIMENSION FX(24),STAGE(24)
      DIMENSION GB(10)
      IF(INDEX6.EQ.1) GO TO 1
      IF(INDEX8.EQ.1) GO TO 2
      IF(INDEX8.LT.4) RETURN
    2 IF(IM.NE.II) RETURN
      DO 15 J=1,24
   15 FX(J)=STAGE(J)
    7 CONTINUE
      DO 4 J=1,7
    4 PDL1(J)=SL1(J)
      PDL1(8)=BET2
      PDL1(9)=AT1
      PDL1(10)=ALFA1*57.296
      PDL1(11)=PPI
      PM=GB(1)
      SGKP=SIGKKP
      GO TO 3
    1 IF(II.NE.IM.OR.ABS(DSIGMA).GT.GB(10)) GO TO 8
      DO 6 J=1,24
    6 FX(J)=STAGE(J)
    9 CONTINUE
      IF(KM.EQ.1) GO TO 7
    3 DO 5 J=1,6
    5 PDL2(J)=SL2(J)
      PDL2(7)=AT1
      PDL2(8)=AT2
      PDL2(9)=ALFA1*57.296
      PDL2(10)=PPI
      PM=GB(1)
      RETURN
    8 IF(II.NE.IN.OR.II.NE.IM) RETURN
      GO TO 9
      END
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
      HQP=HQ
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
      QLA1=QQL(C1A)                                               !压气机入口折合流量
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
     *ALFA2,C2A,AT2,PI,PIK,PIT,LA4,C4A,HT,AKPX)  !各区域参数计算？分离区参数计算
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
      LA2=C2A*UK2/(SIN(ALFA2)*SQRT(2.*9.81*RRR*AK*AT2/(AK+1.)))
      SL1(4)=PI
      SL1(5)=PIK
      SL1(6)=PIT
      SL1(7)=LA2
      IF(INDEX3.NE.0) GO TO 330
      DO 331 J=1,7
  331 PSL1(J)=SL1(J)
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
      CALL BLOCK(II,HQ,QLA1,C1AO,B2O,HTO,ADST,RE,UPR2,
     *RK,HA,RKI,HAI,STAGE,
     *ALFA2,C2A,AT2,PI,PIK,PIT,LA4,C4A,HT,AKPX)
      BETA2=BETA(C2A,R2,ALFA2)
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
      ! 计算出口流量及流速
      CALL AEXIT(C2A,AT2,PI,QLA1,C2AOA,ALFA4,LA4,C4A,II)
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
   49 CALL BLOCK(II,HQ,QLA1,C1AO,B2O,HTO,ADST,RE,UPR2,
     *RK,HA,RKI,HAI,STAGE,
     *ALFA2,C2A,AT2,PI,PIK,PIT,LA4,C4A,HT,AKPX)
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
      DO 302 J=1,7
  302 SL1(J)=PSL1(J)
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
      FUNCTION BETA(C1A,R1,ALFA1)
      IF(C1A.LT.0.01) C1A=0.01
      IF(C1A.GT.0.01) CTGB1=R1/C1A-CTG(ALFA1)
      A=ATAN(CTGB1)
      BETA=1.570796-ATAN(CTGB1)
      RETURN
      END
      FUNCTION FSIGMA(C1A,C1AKP,AL11KP,DDK1)
      XK=10.*(C1A/C1AKP-1.)
      IF(AL11KP.GT.0.8) AL11KP=0.8
      FSIGMA=1.+0.211*(1.-DDK1)*XK*(2*XK+1)/(1.-AL11KP)**2
      RETURN
      END
      SUBROUTINE ALEVEL(HAG,ALFA2,LA2,AK,AMKPA,AM2,K,Z)
      COMMON/A21/KRAN(40)
      REAL HAG,LA2
  801 FORMAT('   ALEVEL   开始')
      IF(KRAN(29).EQ.1) WRITE(6,801)
      AM2=SQRT(2/(AK+1.))*LA2/SQRT(1.-(AK-1.)*LA2*LA2/(AK+1.))
      CALL ELEVEL(HAG,ALFA2,AMKPA)
      Z=ABS(AM2-AMKPA)
      IF(KRAN(35).EQ.1) WRITE(6,*) 'HAG,ALFA2,LA2,AM2,AMKPA,Z'
      IF(KRAN(35).EQ.1) WRITE(6,*) HAG,ALFA2,LA2,AM2,AMKPA,Z
      IF(AM2-AMKPA)1,2,3
    1 K=1
      RETURN
    2 K=2
      RETURN
    3 K=3
      RETURN
      END
      SUBROUTINE ELEVEL(HAG,ALFA2,AMKPA)
      COMMON/A21/KRAN(40)
      REAL HAG
  802 FORMAT('   ELEVEL   HA4AT')
      IF(KRAN(30).EQ.1) WRITE(6,802)
      A=HAG/SIN(ALFA2)
      IF(A-0.825)1,1,2
    1 AMKPA=0.6375*A**2
      RETURN
    2 IF(A-1.09)3,3,4
    3 AMKPA=2.25*A**2-2.62*A+1.064
      RETURN
    4 IF(A-1.3)5,5,6
    5 AMKPA=0.98-2.235*(1.3-A)**2
      RETURN
    6 AMKPA=0.98
      RETURN
      END
      SUBROUTINE ULEVEL(A,AMMAX)
      COMMON/A21/KRAN(40)
  700 FORMAT(6F20.6)
      IF(A.GT.1.081) GO TO 1
      AMMAX=0.78*A*A
      RETURN
    1 IF(A.GT.1.3) GO TO 2
      AMMAX=1.-1.8484*(1.3-A)**2
      RETURN
    2 AMMAX=1.
      RETURN
      END
      SUBROUTINE DOWN1(II,RK,C,C1A,KAG,BETA1,AML1,AMMAX,KM,IM,KK)
      REAL KAG    !叶栅喉道尺寸与步长之比
      DIMENSION RK(20),C(4)
      COMMON/A21/KRAN(40)
      COMMON/A23/INDEX6,INDEX7,INDEX8,INDEX9,GUIDE
      COMMON/A37/WAG1,WAG2,IND8,IND8A
      COMMON/A38/IND8R
  906 FORMAT('   KOLECO ZAPERTO, CTYPEHX  HOMEP',I4)
      IF(RK(16).GT.0) GO TO 1
      AGKA1=KAG/SIN(BETA1)
      CALL ULEVEL(AGKA1,AMMAX)
      IF(AML1.LT.AMMAX) GO TO 2   
    3 IF(KRAN(39).EQ.0) GO TO 10
      WRITE(15,906)II
      WRITE(6,906) II
   10 CONTINUE
      KM=1
      IF(IM.EQ.II) GO TO 21
      IF(IND8R.EQ.0) GO TO 22
      IF(INDEX8.EQ.1) IND8=II
      IF(INDEX8.GE.4) IND8=100*II
   22 CONTINUE
      INDEX8=0
      IF(INDEX6.EQ.0) GUIDE=0
   21 CONTINUE
      IM=II
      KK=2
      RETURN
    1 IF(C1A.GE.C(1)) GO TO 3
    2 KK=1
      RETURN
      END
      SUBROUTINE DOWN2(II,GB,QLA1,LA2,ALFA2,HAG,AK,AML2,AMMAX,KM,IM,KK)
      COMMON/A21/KRAN(40)
      COMMON/A23/INDEX6,INDEX7,INDEX8,INDEX9,GUIDE
      COMMON/A27/SIGMA,DSIGMA,DELTQ
      COMMON/A37/WAG1,WAG2,IND8,IND8A
      COMMON/A38/IND8R
      DIMENSION GB(10)
      REAL LA2
  907 FORMAT('  €弿€悁?噣弲悞, 憭搹厤? 崕寘?',I4)
      AGKA1=HAG/SIN(ALFA2)
      CALL ULEVEL(AGKA1,AMMAX)
      AML2=LA2/SQRT((AK+1.-(AK-1.)*LA2*LA2)/2.)
      IF(AML2.LT.AMMAX) GO TO 1
      IF(KRAN(39).EQ.0) GO TO 20
      WRITE(15,907)II
      WRITE(6,907) II
   20 CONTINUE
      KM=2
      IF(IM.EQ.II) GO TO 21
      IF(IND8R.EQ.0) GO TO 22
      IF(INDEX8.EQ.1) IND8=II
      IF(INDEX8.GE.4) IND8=100*II
   22 CONTINUE
      INDEX8=0
      IF(INDEX6.EQ.0) GUIDE=0
   21 CONTINUE
      IM=II
      KK=2
      RETURN
    1 KK=1
      IF(INDEX6.EQ.1) GO TO 2
      IF(ABS(DELTQ).LE.GB(9)) GO TO 3
    4 IF(GUIDE.GE.0.5) GO TO 5
      RETURN
    3 IF(IM.NE.II) RETURN
      GUIDE=KM+QLA1
      GO TO 4
    5 IF(INDEX8.GE.4) GO TO 6
      IF(INDEX8.NE.0) GO TO 7
      INDEX8=1
    7 CONTINUE
      RETURN
    6 DELTQ=0.25*GB(9)
      RETURN
    2 RETURN
      END
      SUBROUTINE SLOPE1(SL1,II,GKA,C1AKP,BET2,AL11KP,PITKP,PIKKP,F1,
     *C1U,C2A,C2,BETA2,ALFA2,LA2,AT2,ALFA4,LA4,C4A,PI,PIK,PIT,HT,AKPX)
      REAL LA2,LA4
      DIMENSION SL1(7),GKA(4,30)
      COMMON/A2/R1
      COMMON/A6/R2,DK1,DK2,DK4,DDK1,DDK2,DDK4
      COMMON/A7/UK1/A8/ALFA1/A9/RRR,AT1,UK2,UK4/A10/AKH
      COMMON/A13/AK
      COMMON/A15/AKR
      COMMON/A24/IL,IM,IN,KM,KMI
      COMMON/A27R/SIGKKP,SGKP,SMKP
      COMMON/A27/SIGMA,DSIGMA,DELTQ
      COMMON/A21/KRAN(40)
      COMMON/A35/IJ
  404 FORMAT('   SIGMAK SIGMAA SIGMCT PIK PIT PI AKPD')
  700 FORMAT(6F20.6)
  805 FORMAT('   SLOPE1   崁梹')
  905 FORMAT('   SLOPE1   妿崡厤')
      IF(KRAN(34).EQ.1) WRITE(6,805)
      CALL SLOPE3(SL1,GKA,II,C1AKP,BET2,PITKP,PIKKP,AL11KP,F1,
     *HZ,HZ1,HT,PIT,PIK,BETA2,SIGMAK)
      IF(IJ.GT.0) RETURN
      C1A=SL1(1)
      AKR=SQRT(2.*9.81*AK*RRR*AT1/(AK+1.))
      QLA1=QQL(C1A)
      LA2=SL1(7)
      SIGMAO=SL1(4)/SL1(5)
      ELA=(1.-(AK-1.)*LA2*LA2/(AK+1.))**(1./(AK-1.))
      DSITA=(1.-SIGMAO)/(AK*ELA*LA2*LA2/(AK+1.))
       CALL BEXIT(HZ,HZ1,C1A,PIK,QLA1,
     *C1U,C2A,C2U,C2,BETA1,BETA2,ALFA2,LA2,AT2)
      ELA2=(1.-(AK-1.)*LA2*LA2/(AK+1.))**(1./(AK-1.))
      SIGMAA=1.-AK*DSITA*ELA2*LA2*LA2/(AK+1.)
      SIGMCT=SIGMAK*SIGMAA
      PI=PIT*SIGMCT
      FT=AT2-AT1
      DT=ABS(FT)
      IF(DT.LE.0.0001) FT=0.0001*FT/DT
      AKPX=AT1*(PI**((AK-1.)/AK)-1.)/FT
      IF(KRAN(34).EQ.0) RETURN
      WRITE(6,404)
      WRITE(6,700) SIGMAK,SIGMAA,SIGMCT,PIK,PIT,PI,AKPX
      WRITE(6,905)
      RETURN
      END
      SUBROUTINE SLOPE3(SL1,GKA,II,C1AKP,BET2,PITKP,PIKKP,AL11KP,F1,
     *HZ,HZ1,HT,PIT,PIK,BETA2,SIGMAK)
      COMMON/A2/R1
      COMMON/A6/R2,DK1,DK2,DK4,DDK1,DDK2,DDK4
      COMMON/A7/UK1/A8/ALFA1/A9/RRR,AT1,UK2,UK4/A10/AKH
      COMMON/A13/AK
      COMMON/A21/KRAN(40)
      COMMON/A24/IL,IM,IN,KM,KMI
      COMMON/A27/SIGMA,DSIGMA,DELTQ
      COMMON/A27R/SIGKKP,SGKP,SMKP
      COMMON/A35/IJ
      DIMENSION SL1(7),GKA(4,30)
  407 FORMAT('   C1A BETA1 AL11 BET2 PITKP PIKKP AL11KP F,')
  408  FORMAT('   HZ HZ1 HT PIT PIK BETA2 SIGMAK')
  700 FORMAT(6F20.6)
  804 FORMAT('   SLOPE3   崁梹?')
  904 FORMAT('   SLOPE3   妿崡厤')
      C1A=SL1(1)
      BETA1=SL1(2)
      AL11=SL1(3)
      IF(KRAN(32).EQ.0) GO TO 3
      WRITE(6,804)
      WRITE(6,407)
      WRITE(6,700) C1A,BETA1,AL11,BET2,
     *PITKP,PIKKP,AL11KP,F1
      WRITE(6,*) 'C1AKP,SIGMA,ALFA1,AT1,UK1,DDK1,DDK2,DK1,R1,DK2,R2,F1'
      WRITE(6,*) C1AKP,SIGMA,ALFA1,AT1,UK1,DDK1,DDK2,DK1,R1,DK2,R2,F1
    3 IF(IN.NE.II) GO TO 1
      SIGMAK=SIGMA
      GO TO 2
    1 SIGKKP=PIKKP/PITKP
      SIGM=FSIGMA(C1A,C1AKP,AL11KP,DDK1)
      SIGMAK=1.-(1.-SIGKKP)*SIGM
    2 AA=2.*AK*9.81*RRR*AT1/(AK+1.)
      BB=(AK-1.)*R1*(R1-2.*C1A*CTG(ALFA1))*UK1**2/(AK+1.)
      AKP12=AA+BB
      AKP1=SQRT(AKP12)
      AT11=(AK-1.)*(UK1*R1)**2/((AK+1.)*AKP12)
      T11=1+AT11*((DK2*R2/(DK1*R1))**2-1.)
      BETA2=BET2*(1.35-0.35*SIGKKP/SIGMAK)
      QLA11=Q(AL11,AK)
      QQQ=SIGMAK*T11**((AK+1)/(2.*(AK-1.)))
      F2=3.1415*DK2*DK2*(1-DDK2**2)/4.
      QLA21=QLA11*F1*SIN(BETA1)/(F2*SIN(BETA2)*QQQ)
      IF(KRAN(32).EQ.1) WRITE(6,*) 'SIGMAK,SIGKKP,SIGM,AKP1
     *,AKP12,AKP1,AT11,T11,
     *BETA2,QLA11,QQQ,F2,QLA21'
      IF(KRAN(32).EQ.1) WRITE(6,*) SIGMAK,SIGKKP,SIGM,AKP1,AKP12
     *,AKP1,AT11,T11,
     *BETA2,QLA11,QQQ,F2,QLA21
      IF(QLA21.GT.1.AND.QLA21.LT.1.05) QLA21=0.999999
      AL21=RLQMDA(QLA21)
      HH=UK1*R1**2*((DK2*R2/(DK1*R1))**2-1.)
      HHH=R1*SQRT(AKP12)*(AL11*COS(BETA1)-AL21*COS(BETA2)*DK2*
     *R2*SQRT(T11)/(DK1*R1))
      HZ=UK1*(HH+HHH)/9.81
      HZ1=9.81*HZ/(UK1**2)
      HT=HZ1/GKA(3,II)
      PIT=(1.+HZ/(AK*RRR*AT1)*(AK-1.))**(AK/(AK-1.))
      PIK=PIT*SIGMAK
      IF(KRAN(32).EQ.0) RETURN
      WRITE(6,408)
      WRITE(6,700) HZ,HZ1,HT,PIT,PIK,BETA2,SIGMAK
      WRITE(6,904)
      RETURN
      END
      SUBROUTINE SLOPE2(II,SL2,C2AKP,PIKKP,PITKP,PIKP,AL2KP,C2AOA,
     *PI,PIK,PIT,LA4,C4A,ALFA4,AKPX,HT)
      COMMON/A2/R1/A6/R2,DK1,DK2,DK4,DDK1,DDK2,DDK4
      COMMON/A7/UK1/A8/ALFA1
      COMMON/A9/RRR,AT1,UK2,UK4/A10/AKH
      COMMON/A11/C2AO,C4AO,ALFA4O/A13/AK
      COMMON/A21/KRAN(40)
      COMMON/A24/IL,IM,IN,KM,KMI
      COMMON/A27/SIGMA,DSIGMA,DELTQ
      COMMON/A29/AT2
      DIMENSION SL2(6)
      REAL LA4
  405 FORMAT('   QLA1 C2A AT1 AT2 PIK PIT')
  406 FORMAT('   PI SIGMAK SIGMAA SIGMCT ALFA4 C4A')
  803 FORMAT('   SLOPE2   崁梹?')
  903 FORMAT('   SLOPE2   妿崡厤')
  700 FORMAT(6F20.6)
      IF(KRAN(31).EQ.1) WRITE(6,803)
      QLA1=SL2(1)
      C2A=SL2(2)
      AT2=SL2(3)
      PIK=SL2(4)
      PIT=SL2(5)
      PI=SL2(6)
      IF(KRAN(31).EQ.1) WRITE(6,405)
      IF(KRAN(31).EQ.1) WRITE(6,700) QLA1,C2A,AT1,AT2,PIK,PIT
      IF(IN.EQ.II) GO TO 1
    3 CONTINUE
      SIGAKP=PIKP/PIKKP
      SIGM=FSIGMA(C2A,C2AKP,AL2KP,DDK1)
      SIGMAA=1.-(1.-SIGAKP)*SIGM
      GO TO 2
    1 IF(KMI.NE.2) GOTO 3
      SIGMAA=SIGMA
    2 CONTINUE
      SIGMAK=PIK/PIT
      SIGMCT=SIGMAK*SIGMAA
      PI=PIT*SIGMCT
      FT=AT2-AT1
      DT=ABS(FT)
      IF(DT.LE.0.0001) FT=0.0001*FT/DT
 1111 FORMAT(10F12.6)
      AKPX=AT1*(PI**((AK-1.)/AK)-1.)/FT
      ! 计算出口流量及流速
      CALL AEXIT(C2A,AT2,PI,QLA1,C2AOA,ALFA4,LA4,C4A,II)
      IF(KRAN(31).EQ.0) RETURN
      WRITE(6,406)
      WRITE(6,700) PI,SIGMAK,SIGMAA,SIGMCT,ALFA4,C4A
      WRITE(6,903)
      RETURN
      END
      SUBROUTINE BEXIT(HZ,HZ1,C1A,PIK,QLA1,
     *C1U,C2A,C2U,C2,BETA1,BETA2,ALFA2,LA2,AT2)
      REAL LA2
      COMMON/A2/R1
      COMMON/A6/R2,DK1,DK2,DK4,DDK1,DDK2,DDK4
      COMMON/A7/UK1/A8/ALFA1/A9/RRR,AT1,UK2,UK4/A10/AKH
      COMMON/A13/AK
      COMMON/A21/KRAN(40)
  409 FORMAT('   HZ HZ1 C1A PIK QLA1 AT1')
  410 FORMAT(' = AT2 BETA1 C1U C2A C2 ALFA2 BETA2 LA2')
  700 FORMAT(6F20.6)
  505 FORMAT('   EXIT   开始')
  605 FORMAT('   EXIT   结束')
      IF(KRAN(33).EQ.1) WRITE(6,505)
      IF(KRAN(33).EQ.1) WRITE(6,409)
      IF(KRAN(33).EQ.1) WRITE(6,700) HZ,HZ1,C1A,PIK,QLA1,AT1
      AT2=AT1+HZ*(AK-1.)/(AK*RRR)
      C1U=C1A*CTG(ALFA1)
      C2U=(HZ1+R1*C1U)/(R2*((DK2)/(DK1))**2)
      T2A=AT2-(C2U*UK2)**2*(AK-1.)/(2.*AK*RRR*9.81)
      AKR2A=SQRT(2*AK*9.81* RRR*T2A/(AK+1.))
      QQ1=QLA1*SIN(ALFA1)/PIK
      QQ2=(DK1/DK2)**2*((1.-DDK1**2)/(1.-DDK2**2))
      QQ3=SQRT(AT2/AT1)*(AT2/T2A)**((AK+1.)/(2.*(AK-1.)))
      QLA2A=QQ1*QQ2*QQ3
      ALA2A=RLQMDA(QLA2A)
      C2A=ALA2A*AKR2A
      C2=SQRT(C2A**2+(C2U*UK2)**2)
      C2A=C2A/UK2
      TAN2=C2A/C2U
      ALFA2=1.570796-ARCCTG(TAN2)
      BETA1=BETA(C1A,R1,ALFA1)
      BETA2=BETA(C2A,R2,ALFA2)
      AKR2=SQRT(2.*9.81*AK*RRR*AT2/(AK+1.))
      LA2=C2/AKR2
      IF(KRAN(33).EQ.0) RETURN
      WRITE(6,410)
      WRITE(6,700) AT2,BETA1,C1U,C2A,C2,ALFA2,BETA2,LA2
      WRITE(6,605)
      RETURN
      END
      SUBROUTINE PARINV(X,A,F,N,R)    !对应于不同流道形式
      DIMENSION A(N),F(N)
        IF(X.LT.A(1)) GOTO 11
        IF(X.GT.A(N)) GOTO 4
        K1=1
        K2=N
    2   K3=K2-K1
        IF(K3.LE.1) GOTO 6    ! N小于3 GOTO 6
        K3=K1+K3/2
        IF(A(K3)-X) 7,8,9
    7   K1=K3
        GOTO 2
    9   K2=K3
        GOTO 2
    8   R=F(K3)
        RETURN
    3   B1=A(K1)
        B2=A(K1+1)
        B3=A(K1+2)
        B4=F(K1)
        B5=F(K1+1)
        B6=F(K1+2)
        R=B4*((X-B2)*(X-B3))/((B1-B2)*(B1-B3))+B5*((X-B1)*(X-B3))
     *  /((B2-B1)*(B2-B3))+B6*((X-B1)*(X-B2))/((B3-B1)
     *  *(B3-B2))
        RETURN
    6   IF(K2.NE.N) GOTO 3
        K1=N-2
        GOTO 3
    4   C=ABS(X-A(N))
        IF(C.LT.1.E-7) GOTO 5
        K1=N-2
   13   WRITE(6,41) X
   41   FORMAT(25H X IS OUT OF THE INTERVAL,3H X=,F15.9)
        R=0.
        RETURN
    5   R=F(N)
        RETURN
   11   C=ABS(X-A(1))
        IF(C.LT.0.1E-7) GOTO 12
        K1=1
        GOTO 13
   12   R=F(1)
        RETURN
      END
      ! 计算出口流量及流速
      SUBROUTINE AEXIT(C2A,AT2,PI,QLA1,
     *C2AOA,ALFA4,LA4,C4A,II)  !静叶出口计算
      REAL LA4,NUA
      COMMON/A2/R1/A6/R2,DK1,DK2,DK4,DDK1,DDK2,DDK4
      COMMON/A7/UK1/A8/ALFA1
      COMMON/A9/RRR,AT1,UK2,UK4/A10/AKH
      COMMON/A11/C2AO,C4AO,ALFA4O/A13/AK
      COMMON/A42/ALF4C
      COMMON/BIS/BIS(30)
      AKR2=SQRT(2.*9.81*RRR*AK*AT2/(AK+1.))
      WW1=DK1**2/DK4**2
      WW2=(1-DDK1**2)/(1-DDK4**2)
      WW3=SQRT(AT2/AT1)
      WW4=SIN(ALFA1)/SIN(ALFA4O)
      WW5=1.-BIS(II)
      QLA4=QLA1*WW1*WW2*WW3*WW4/PI    ! 出口折合流量计算
      LA4=RLQMDA(QLA4)
      C4A=LA4*AKR2*SIN(ALFA4O)/UK4
      ! 满足以下条件之后
      IF(ALFA4O.GE.1.57096.OR.ALF4C.GT.0.) THEN
          ALFA4=ALFA4O
          RETURN
      END IF
      IF(C2A.LT.0.01) C2A=0.01
      NUA=C2A/C2AOA
      IF(NUA.GE.1.0) THEN
          DHA=0.1*(1-NUA)
      ELSE
          DHA=0.15*(1-NUA)
      END IF
      WW1=R2/C2A
      WW2=C4AO/C2AO
      IF(C4A.LT.0.01) C4A=0.01
      WW3=C2A/C4A
      WW4=COS(ALFA4O)/SIN(ALFA4O)
      CTA4=WW3*(WW1*DHA+WW2*WW4)  !静叶出口气流角余切
      ALFA4=1.570796-ATAN(CTA4)
      QLA4=QLA4*SIN(ALFA4O)/SIN(ALFA4)
      LA4=RLQMDA(QLA4)
      C4A=LA4*AKR2*SIN(ALFA4)/UK4
      RETURN
      END
      ! 离心压气机初步计算
      SUBROUTINE CENTRO(BETIN,ALF1CIN,DBT1IN,DP1IN,NIN,GIN,T0ZIN,P0ZIN,
     *PICTIN,IERR,D2IN,H2IN,LPKIN,DGABIN,P6ZIN,T6ZIN,KPDCBKIN,
     *H3IN,D3IN,H4IN,D4IN,D5IN,LAM3IN,PSI,GAMLD)
      IMPLICIT REAL (K-N)
      IMPLICIT INTEGER (Z)
      INTEGER ANSWER
      PARAMETER(PJ=3.14129,DEG=57.2958)
      OPEN(8,FILE='CENTRO.OUT')
      WRITE(8,995)
 995  FORMAT(1X,'THE CENTRIFUGAL COMPRESSOR PRELIMINARY CALCULATION')
      WRITE(8,994)
 994  FORMAT(1X,'COPYRIGHT (C) DMITRY M. AFANASSIEV , CIAM ,  1992 ')
      WRITE(8,112)
      WRITE(8,991) BETIN*DEG,ALF1CIN*DEG
 991  FORMAT(1X,'BETTA 2L=',F6.2,5X,'ALPHA 1=',F6.2)
      WRITE(8,992) DBT1IN,DP1IN
 992  FORMAT(1X,'HUB DIAM.=',F6.3,4X,'TIP DIAM.=',F7.3)
      WRITE(8,993) NIN,GIN,T0ZIN,P0ZIN,PICTIN
 993  FORMAT(1X,'N=',F8.1,3X,'G=',F6.3,2X,'T0 TOTAL=',F7.2,2X,
     *'P0 TOTAL=',F9.1,4X,/,1X,'PI NEED=',F6.3)
      IERR=0
      AL1=ALF1CIN*DEG
      N=NIN
      G=GIN
      BETTL2=BETIN
      ALF1C=ALF1CIN
      ALF1P=ALF1CIN
      ALF1B=ALF1CIN
      DBT1=DBT1IN
      DP1=DP1IN
      T0Z=T0ZIN
      P0Z=P0ZIN
      PICT=PICTIN
      IF(PICT.LE.1.1) THEN
      WRITE(8,*) 'NO CENTRIFUGAL STAGE: TOO LOW PRESSURE RATIO=',PICT
      IERR=1
      END IF
      IF(PICT.GT.20.) THEN
      WRITE(8,*) 'NO CENTRIFUGAL STAGE: TOO GREAT PRESSURE RATIO=',PICT
      IERR=1
      END IF
      IF(AL1.LE.25.) THEN
      IERR=1
      WRITE(8,*) 'NO CENTRIFUGAL STAGE: TOO LOW ALPHA1=',AL1
      END IF
      DBTODPR = DBT1/DP1
      IF(DBTODPR.LT.0.2) THEN
      WRITE(8,*) 'NO CENTRIFUGAL STAGE: TOO LITTLE HUB RATIO=',DBTODPR
      IERR=1
      END IF
      IF(DBTODPR.GT.0.9) THEN
      WRITE(8,*) 'NO CENTRIFUGAL STAGE: TOO GREAT HUB RATIO=',DBTODPR
      IERR=1
      END IF
      IF(IERR.EQ.1) GOTO 1
      F1=(DP1*DP1-DBT1*DBT1)*PJ/4.
      DC1=SQRT((DP1*DP1+DBT1*DBT1)/2.)
      OMEGA=PJ*N/30.
      UC1=OMEGA*DC1/2.
      IF(UC1.LT.10.) THEN
      WRITE(8,*) 'NO CENTRIFUGAL STAGE: TOO LITTLE U1=',UC1
      IERR=2
      END IF
      K=1.4
      R=287.14
      KG1=1.03
      MCR=SQRT(K*(2/(K+1))**((K+1)/(K-1)))
      M=MCR/SQRT(R)
      QOTL1C=G*KG1*SQRT(T0Z)/(M*SIN(ALF1C)*P0Z*F1)
      IF(QOTL1C.LT.0.0001) THEN
      WRITE(8,*) 'NO CENTRIFUGAL STAGE:  NO INLET VELOCITY'
      IERR=2
      GOTO 1
      END IF
      IF(QOTL1C.GT.1.0001) THEN
      WRITE(8,*) 'NO CENTRIFUGAL STAGE: TOO LITTLE INLET AREA'
      WRITE(8,*) 'Q(K,LAM1)=',QOTL1C
      IERR=2
      GOTO 1
      END IF
      IF(NS.LT.10.OR.NS.GT.70.) THEN
      WRITE(8,*) 'NO CENTRIFUGAL STAGE: NS=',NS
      IERR=2
      GOTO 1
      END IF
      CAPC1=1.03
      CABC1=1.03
      BETTAF=1.25
      KG2=1.05
      KD2=1.04
      LATP=0.022
      LATP1=0.01
      DEL4=0.
      ZLOP=0
      KG3=1.05
      KG4=1.05
      KG5=1.03
      ALFA6=PJ/2.
      SIG34=0.9
      SIG45=0.995
      SIG56=0.995
      KG6=1.05
      SIG01=1.
      IF(PICT.LE.8.) THEN
      U2=-4.441*PICT*PICT+93.73*PICT+172.372
      ELSE
      U2=24.8*PICT+440.
      END IF
      U2=U2*1.2
      D2=60.*U2/PJ/N
      KPDAD=0.920-(30.*PICT+0.1093*(BETTL2*DEG)**2-
     *8.944*(BETTL2*DEG))/10000.
      IF(DBTODPR.GT.0.6)KPDAD=KPDAD-(10*DBTODPR-6)/100
      GPRIVED=G*101325./P0Z*SQRT(T0Z/288.)
      IF(GPRIVED.LT.0.7) KPDAD=KPDAD*(0.253*GPRIVED+0.807)
      IF(GPRIVED.GE.0.7.AND.GPRIVED.LT.2.)
     *KPDAD=KPDAD*(0.0138*GPRIVED+0.972)
      IF(NS.LT.35)KPDAD=KPDAD+0.001466*(NS-35)
      IF(NS.GT.45)KPDAD=KPDAD-(NS-45.)/5000
      LA1C=0.
      BETT1B=PJ/2.
      KPDLA=KPDAD
      PR2=0.467-DBT1/DP1/3
      ZSCORE=1
      PR=1.
      PPR=1.
 100  CONTINUE
      IF(LA1P.GT.0.5) KPDAD=KPDLA+3*(1-2*LA1C)/100
      H2=D2*(0.0016*NS-0.021)
      H34=H2
      D3=D2*1.1
      D4=D2*1.38
      D5H=D2*1.55
      F34=2.5
      H4=H34+D2*0.0122468
      D6H=D5H
      D6BH=D6H-1.8*H2
      D5BH=D6BH
      L45=19./282.*D2
      ZPKMAX=PJ*DBT1*SIN(BETT1B)/0.003
      IF(ZPKMAX.LT.ZPK) ZPK=ZPKMAX
      MI=0.985*(1.4707+0.01717*BETTL2*DEG)/(3.333-1./TAN(BETTL2))
      MI=MI*(0.9253+0.00267*ZPK)
      IF(MI.GT.0.97) MI=0.97
      IF(MI.LT.0.6) MI=0.6
      GOTO 801
 802  CONTINUE
      IF(ANSWER.EQ.3) THEN
      WRITE(8,*) 'NS=',NS
      WRITE(8,*) 'NO CENTRIFUGAL STAGE: IMPOSSIBLE PARAMETERS'
      IERR=3
      GOTO 1
      END IF
      SIG34=1.098-(0.16+0.006*F34)*(LAM3**0.85)
      SIG45=-0.03125*LAMDA4+1.001875
      SIG56=-0.0385*LAMDA5+1.0034
      DEVIA=PIKAZ-PICT
      ODEVIA=ABS(DEVIA)/PIKAZ
      WRITE(8,*)ZSCORE,'     PICT=',PIKAZ,'       ',ODEVIA
      IF(DEVIA.LT.-0.005.AND.ZSCORE.LE.90.AND.PIKAZ.NE.PPR) THEN
      PPR=PR
      PR=PIKAZ
      ZSCORE=ZSCORE+1
      D2=D2*(ODEVIA/10.+1.)
      GOTO 100
      END IF
      IF(DEVIA.GT.0.005.AND.ZSCORE.LE.90.AND.PIKAZ.NE.PPR) THEN
      PPR=PR
      PR=PIKAZ
      ZSCORE=ZSCORE+1
      D2=D2/(ODEVIA/10.+1.)
      GOTO 100
      END IF
      IF(ZSCORE.GT.90) THEN
      WRITE(8,*)'NO CENTRIFUGAL STAGE: TOO HARD PARAMETERS'
      GOTO 1
      END IF
      WRITE(8,*)'WARNING: THIS VARIANT IS PRELIMINARY ONE'
      WRITE(8,101)G,N,R,T0Z,CAPC1
      WRITE(8,102)GPRIVED,NPRIVED,K,P0Z,CABC1
      WRITE(8,112)
      WRITE(8,103)DP1,DBT1,D5H,D5BH,D6H
      WRITE(8,104)DVTOD1,D1OD2,D3OD2,D4OD3,D6BH
      WRITE(8,105)GAMBLD*DEG,GAMLD*DEG,F34,H2OD2,D5HOD2
      WRITE(8,112)
      WRITE(8,106)KD2,BETTL2*DEG,BETTAF,ALFAF,KPDAD
      WRITE(8,107)DZIT12,LATP,LATP1,MI,OMEGA
      WRITE(8,113)DELOT*DEG,L45,POWER
      WRITE(8,112)
      WRITE(8,108)      NS,PIKAZ,KPDADK
      WRITE(8,109)SIGBC,C1ACOTN,NSND,HZ,HADKZ1
      WRITE(8,112)
      WRITE(8,114)REACT,REACTD,ZPK
101   FORMAT(3X,'G=',F6.3,5X,'N=',F8.1,3X,'R=',F7.3,1X,'T0*=',F6.2,
     *3X,'C1APER/CP=',F4.2)
102   FORMAT(1X,'GPR=',F6.3,3X,'NPR=',F8.1,3X,'K=',F6.4,2X,'P0*=',F8.1,
     *2X,'C1ACP/BT=',F4.2)
103   FORMAT(2X,'D1=',F5.4,5X,'DBT=',F5.4,2X,'D5HAP=',F5.3,3X,'D5BH=',
     *F5.3,3X,'D6HAP=',F5.3)
104   FORMAT(1X,'DBT/D1=',F4.3,2X,'D1/D2=',F4.3,2X,'D3/D2=',F5.3,3X,
     *'D4/D3=',F5.3,2X,'D6BH=',F5.3)
105   FORMAT(1X,'GAMBLD=',F5.2,1X,'GAMLD=',F5.2,3X,'F34=',F4.2,4X,
     *'H2/D2=',F4.3,3X,'GABARIT=',F5.3)
106   FORMAT(5X,'KD2=',F4.2,2X,'BETT2L=',F5.1,2X,'BETTF=',F4.2,2X,
     *'ALFAF=',F4.3,1X,'KPD2=',F4.3)
107   FORMAT(2X,'DZIT12=',F5.3,1X,'LAMTP=',F4.3,3X,'LAMTP1=',F4.3,2X,
     *'MI=',F5.4,2X,'OMEGA=',F8.2)
108   FORMAT(3X,'      ',6X,'NS=',F5.2,3X,
     *'PICT=',F5.2,4X,'KPDAD=',F4.3,2X)
109   FORMAT(3X,'SIGBC=',F4.3,2X,'KRAS=',F5.3,1X,'NSND=',F4.2,5X,
     *'HZ=',F4.2,2X,'HADCT=',F4.2)
112   FORMAT(3X)
113   FORMAT(3X,'DELTOTC=',F5.2,12X,'L45=',F4.3,5X,'POWER=',F11.2)
114   FORMAT(3X,'ROK=',F5.3,4X,'ROKD=',F5.3,14X,'Z=',I3)
201   FORMAT(6X,'         ',3X,'   TIP   ',4X,'   MID    ',4X,'  HUB '
     *,6X,' EXIT')
202   FORMAT(6X,' D M  ',4(6X,F7.4))
203   FORMAT(6X,'U  M/C',4(6X,F7.3))
204   FORMAT(6X,'P*  KP',4(6X,F7.1))
205   FORMAT(6X,'LAMDA ',4(6X,F7.4))
206   FORMAT(6X,'LAMDM ',4(6X,F7.4))
207   FORMAT(6X,'LAMDU ',4(6X,F7.4))
208   FORMAT(6X,'C  M/C',4(6X,F7.3))
209   FORMAT(6X,'CM  M/C',5X,F7.3,6X,F7.3,6X,F7.3,6X,F7.3)
210   FORMAT(6X,'CU  M/C',5X,F7.3,6X,F7.3,6X,F7.3,6X,F7.3)
211   FORMAT(6X,'W   M/C',5X,F7.3,6X,F7.3,6X,F7.3,6X,F7.3)
212   FORMAT(4X,'BET  GRAD',6X,F5.2,8X,F5.2,8X,F5.2,8X,F5.2)
213   FORMAT(6X,'LAM  OTH',5X,F5.3,8X,F5.3,8X,F5.3,8X,F5.3)
214   FORMAT(7X,'P  KP.',3X,F8.1,5X,F8.1,5X,F8.1,5X,F8.1)
215   FORMAT(6X,'T   K',8X,F5.1,8X,F5.1,8X,F5.1,8X,F5.1)
216   FORMAT(6X,'ALFA ',7X,F6.2,7X,F6.2,7X,F6.2,7X,F6.2)
218   FORMAT(5X,'W1/W2 ',7X,F6.4,7X,F6.4,7X,F6.4)
      WRITE(8,201)
      WRITE(8,202)DP1,DC1,DBT1,D2
      WRITE(8,203)UP1,UC1,UB1,U2
      WRITE(8,204)P1Z/1000,P1Z/1000,P1Z/1000,P2Z/1000
      WRITE(8,205)LAM1P,LAM1C,LAM1B,LAMDA2
      WRITE(8,206)LAM1AP,LAM1AC,LAM1AB,LAM2M
      WRITE(8,207)LAM1UP,LAM1UC,LAM1UB,LAM2U
      WRITE(8,208)C1P,C1C,C1B,C2
      WRITE(8,209)C1AP,C1AC,C1AB,C2R
      WRITE(8,210)C1UP,C1UC,C1UB,CU2
      WRITE(8,211)W1P,W1C,W1B,W2
      WRITE(8,212)BETT1P*DEG,BETT1C*DEG,BETT1B*DEG
     *,BETTA2*DEG
      WRITE(8,213)LA1P,LA1C,LA1B,LAM2W
      WRITE(8,214)P1P/1000,P1C/1000,P1B/1000,P2/1000
      WRITE(8,215)T1P,T1C,T1B,T2
      WRITE(8,216)ALF1P*DEG,ALF1C*DEG,ALF1B*DEG
     *,ALFA2*DEG
      WRITE(8,218)W1P/W2,W1COW2,W1B/W2
355   FORMAT(19X,'1',8X,'2',8X,'3',8X,'4',8X,'5',8X,'6')
356   FORMAT(9X,'H M',3X,6(F6.4,3X))
357   FORMAT(9X,'D M',3X,6(F6.4,3X))
3571  FORMAT(9X,'KG ',3X,6(F5.3,4X))
358   FORMAT(9X,'T* K',3X,6(F6.1,3X))
359   FORMAT(7X,'P* KP',2X,6(F8.1,1X))
360   FORMAT(8X,'T  K',3X,6(F6.1,3X))
361   FORMAT(6X,'P  KP',3X,6(F8.1,1X))
362   FORMAT(5X,'RO KG/M3',3X,6(F5.2,4X))
363   FORMAT(5X,'U M/CEK',4X,2(F5.1,4X))
364   FORMAT(7X,'CM M/CEK',1X,6(F6.2,3X))
365   FORMAT(7X,'CU M/CEK',1X,6(F6.2,3X))
565   FORMAT(7X,'C  M/CEK',1X,6(F6.2,3X))
366   FORMAT(7X,'LAMDA   ',1X,6(F6.4,3X))
367   FORMAT(7X,'LAMDAM  ',1X,6(F6.4,3X))
368   FORMAT(7X,'LAMDAU  ',1X,6(F6.4,3X))
384   FORMAT(7X,'ALFA GR  ',6(F6.2,3X))
383   FORMAT(7X,'PI *    ',9X,5(F6.3,3X))
382   FORMAT(7X,'SIGMA   ',1X,9X,5(F5.3,4X))
371   FORMAT(7X,'HAD J/KG ',7X,6(F8.0,1X))
      WRITE(8,355)
      WRITE(8,356)H1,H2,H34,H4,H5,H6
      WRITE(8,357)DC1,D2,D3,D4,DC5,D6C
      WRITE(8,3571)KG1,KG2,KG3,KG4,KG5,KG6
      WRITE(8,358)T0Z,T2Z,T2Z,T2Z,T2Z,T2Z
      WRITE(8,359)P1Z/1E3,P2Z/1E3,P3Z/1E3,P4Z/1E3,P5Z/1E3,P6Z/1E3
      WRITE(8,360)T1C,T2,T3,T4,T5,T6
      WRITE(8,361)P1C/1E3,P2/1E3,P3/1E3,P4/1E3,P5/1E3,P6/1E3
      WRITE(8,362)RO1,RO2,RO3,RO4,RO5,RO6
      WRITE(8,364)C1AC,C2R,CM3,CR4,CA5,CA6
      WRITE(8,365)C1UC,CU2,CU3,CU4,CU5,CU6
      WRITE(8,565)C1C,C2,C3,C4,C5,C6
      WRITE(8,366)LAM1C,LAMDA2,LAM3,LAMDA4,LAMDA5,LAMDA6
      WRITE(8,367)LAM1AC,LAM2M,LAM3M,LAMR4,LAA5,LAA6
      WRITE(8,368)LAM1UC,LAM2U,LAM3U,LAMU4,LAMU5,LAMU6
      WRITE(8,384)ALF1C*DEG,ALFA2*DEG,ALFA3*DEG,
     *ALFA4*DEG,ALFA5*DEG,ALFA6*DEG
      WRITE(8,383)PI2Z,PI3Z,PI4Z,PI5Z,PIKAZ
      WRITE(8,382)      SIG12,SIG23,SIG34,SIG45,SIG56
      WRITE(8,371)HAD2,HAD3,HAD4,HAD5,HADKZ
  1   CLOSE(8)
      D2IN=D2
      H2IN=H2
      LPKIN=PR2*D2
      DGABIN=D5H
      P6ZIN=P6Z
      T6ZIN=T2Z
      KPDCBKIN=KPDADK
      H3IN=H34
      D3IN=D3
      H4IN=H4
      D4IN=D4
      D5IN=D5BH
      LAM3IN=LAM3
      RETURN
 801  CONTINUE
      ANSWER=0
      CP=K*R/(K-1)
      GPRIVED=G*101330*SQRT(T0Z/288.14)/P0Z
      NPRIVED=N*SQRT(288./T0Z)
      DVTOD1=DBT1/DP1
      D1OD2=DP1/D2
      H2OD2=H2/D2
      D3OD2=D3/D2
      D4OD3=D4/D3
      GAMBLD=ATAN(2*(H34-H2)/(D3-D2))
      GAMLD=ATAN(2*(H4-H34)/(D4-D3))
      D5HOD2=D5H/D2
      D6HOD2=D6H/D2
      H1=(DP1-DBT1)/2.
      OMEGA=PJ*N/30.
      F1=(DP1*DP1-DBT1*DBT1)*PJ/4.
      DC1=SQRT((DP1*DP1+DBT1*DBT1)/2.)
      UC1=OMEGA*DC1/2.
      UP1=OMEGA*DP1/2.
      UB1=OMEGA*DBT1/2.
      P1Z=P0Z*SIG01
      MCR=SQRT(K*(2/(K+1))**((K+1)/(K-1)))
      M=MCR/SQRT(R)
      QOTL1C=G*KG1*SQRT(T0Z)/(M*SIN(ALF1C)*P1Z*F1)
      IF(ABS(ALF1P-PJ/2).LT.1E-3)THEN
      C1UP=0.
      ELSE
      C1UP=C1AP/TAN(ALF1P)
      END IF
      IF(ABS(ALF1B-PJ/2).LT.1E-3)THEN
      C1UB=0.
      ELSE
      C1UB=C1AB/TAN(ALF1B)
      ENDIF
      C1P=SQRT(C1AP*C1AP+C1UP*C1UP)
      C1B=SQRT(C1AB*C1AB+C1UB*C1UB)
      LAM1P=C1P/ACR1
      LAM1B=C1B/ACR1
      LAM1AP=C1AP/ACR1
      LAM1AB=C1AB/ACR1
      LAM1UP=C1UP/ACR1
      LAM1UB=C1UB/ACR1
      W1UP=UP1-C1UP
      W1UC=UC1-C1UC
      W1UB=UB1-C1UB
      BETT1P=ATAN(C1AP/W1UP)
      BETT1C=ATAN(C1AC/W1UC)
      BETT1B=ATAN(C1AB/W1UB)
      W1P=W1UP/COS(BETT1P)
      W1C=W1UC/COS(BETT1C)
      W1B=W1UB/COS(BETT1B)
      RO1=P1C/R/T1C
      U2=OMEGA*D2/2.
      F2=PJ*D2*H2
      A1=C1UC*UC1/(U2*U2)
      ALFAF=0.04
      A0=0.35
      A21=TAN(BETTL2)
13    A2=1.-A0/A21
      HAD21=KPDAD*(MI*A2+ALFAF-A1)
      HAD2=HAD21*U2*U2
      IF(HAD2.LE.0.001) THEN
      WRITE(8,*) ' TOO LITTLE WORK  FACTOR'
      END IF
      END
