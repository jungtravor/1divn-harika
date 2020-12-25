﻿
      SUBROUTINE BLOCK(II,HQ,QLA1,C1AO,B2O,HTO,ADST,RE,UPR2,RK,
     *HA,RKI,HAI,STAGE,ALFA2,C2A,AT2,PI,PIK,
     *PIT,LA4,C4A,HT,AKPX)
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
      CALL YPPKC(HQ,C1AO,HTO,B2O,AML1,ADST,PIK,ALFA2,LA2,C2A,
     *C2U,AT2,II,RK,HA,RKI,HAI,STAGE,AKP,XFK,XFA,
     *UPR2,PI,PIT,ALFA4,LA4,C4A,HT,AKPX)
      IF(HQC.LE.0.) THEN
          IF(RK(20).GT.1E-6) THEN         ! 输入参数中提供：RK(20)分离边界上的空气动力载荷判据值
              FGK=RK(20)
              FGKA=FGK*1.2
          ELSE
              CALL FGCA(ANT,RE,FGK,FGKA)  ! 计算分离判据
          END IF
          IF(INDEX5.GT.0) RETURN
          XXX=AK*RRR*9.81
          HQCK=HQCP(C1AO,B2O,XXX,FGK)     ! 分离边界点相对流量系数的初值确定
          CALL OGCC(HQ,HQCK,C1AO,HTO,B2O,AML1,ADST,FGKA,
     *    II,RK,HA,RKI,HAI,STAGE,
     *    AKP,XFK,XFA,UPR2,HQC,ALIN)
      END IF
      IF(HQC.GT.50.) RETURN
      IF(KRAN(39).EQ.1) WRITE(7,9992)II,FGK,FGKA
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
      IF(KRAN(6).NE.0) WRITE(7,933) HC,H,HT,VTAC,VTAT,AKPX
      CALL PPZPKC(HT,AKPX,QLA1,C1A,C2AOA,PIK,ALFA2,LA2,C2A,C2U,AT2,
     *PI,PIT,ALFA4,LA4,C4A,II)
      IF(KRAN(6).NE.0) WRITE(7,933) HQ,C1A,PIK,ALFA2,LA2,C2A,
     *C2U,AT2,PI,PIT,ALFA4,LA4,C4A,HT,AKPX
      END