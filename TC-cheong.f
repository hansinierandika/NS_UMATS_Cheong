C
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,STRAN,DSTRAN,
     2 TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,MATERL,NDI,NSHR,NTENS,
     3 NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,CELENT,
     4 DFGRD0,DFGRD1,NOEL,NPT,KSLAY,KSPT,KSTEP,KINC)
C
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*8 MATERL	
C      
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),
     4 DFGRD0(3,3),DFGRD1(3,3)
C
C
      DIMENSION A(6),STRI1(6),AF(6)
      DIMENSION ESTRS(6),DE(6,6),DSUM(6,6)
C
C***********************************************************************
C
C UMAT FOR Mike Jefferies's NORSAND, GEOTECHNIQUE(1993)
C UMAT WAS CODED BY GANESH DASARI, CAMBRIDGE UNIVERSITY
C BEEN MODIFIED BY TZI PIAU CHEONG, CAMBRIDGE UNIVERSITY - NOVEMBER 2004
C 1. MATSUOKA & NAKAI'S FAILURE CRITERIA
C 2. Emax AND Emin CALCULATION
C 3. EXPONENTIAL HARDNEING MODULUS AND INTRODUCED OF M/Mtc COEFFICIENT
C 4. HARDNEING RULE, DIFFERENTIAL w.p.t. DEVIATORIC STRAIN 
C
C***********************************************************************
C
C	FOLLOWING LINES SHOW HOW TO USE THIS MODEL IN ABAQUS INPUT FILE
C
C	*SOLID SECTION,ELSET=ALLE,MATERIAL=NORS
C	*MATERIAL,NAME=NORS
C	*USER MATERIAL,CONSTANTS=11
C        RM,EMAX,EMIN,RNOV,HARD,GS,EXP,POIS
C        DILC,ISWIT,TOL
C
C	 RM      Critical Stress Ratio in tx compression, M
C	 EMAX    Maximum void ratio, emax
C      EMIN    Minimum void ratio, emin
C      RNOV    Value of N in Nova'f flow rule, N
C      HARD    Hardening Modulus, h
C      GS      Shear modulus multiplier, A
C	 EXP     Pressure Exponnent to calculate GS, n		
C      POIS    Poisson's Ratio, v
C      DILC    Maximum Dilation Coefficient, about 3.5
C      ISWIT   1.0 constant H, any other number for exponential H
C      TOL     TOLERENCE FOR STRESS INTEGRATION (1E-2 TO 1E-5)
C        
C STATEV(1)=  VOID RATIO
C STATEV(2)=  PMAX STRESS
C STATEV(3)=  STRESS STATE POISITION WRT YSF 0,1,2,3,4
C STATEV(4)=  VALUE OF Q
C STATEV(5)=  STRESS RATIO Q/P
C STATEV(6)=  DEVIATORIC STRAIN
C STATEV(7)=  VOLUMETRIC STRAIN
C STATEV(8)=  VALUE OF PSI
C STATEV(9)=  PLASTIC COMPONENT OF DEVIATORIC STRAIN
C STATEV(10)= PLASTIC COMPONENT OF VOLUMETRIC STRAIN
C STATEV(11)= VALUE OF PI
C STATEV(12)= VALUE OF PI_MAX
C STATEV(13)= VALUE OF THETA
C STATEV(14)= VALUE OF M (BASED ON LODE ANGLE,THETA)
C STATEV(15)= -1.0 IF EIGEN VALUE IS NEGATIVE, OTHERWISE 0.0
C STATEV(16)= 0.0 ELASTIC, 1.0 PLASTIC
C
C INITIALISE LOCAL VARIABLES
C
C
	CALL ZERO1(A,6)
	CALL ZERO1(AF,6)
	CALL ZERO1(ESTRS,6)
	CALL ZERO2(DE,6)
	CALL ZERO2(DSUM,6)
C
C	
	NDIM=3
	IF(NTENS.EQ.4)NDIM=2
C
C     DEFINE MATERIAL PARAMETERS
C
      RM=PROPS(1)
      EMAX=PROPS(2)
      EMIN=PROPS(3)
      RNOV=PROPS(4)
      HARD=PROPS(5)
	RA=PROPS(6)
	REXPN=PROPS(7)
	RV=PROPS(8)
      DILC=PROPS(9)
      ISWIT=PROPS(10)
      TOL=PROPS(11)
C
C     INITIAL PMAX AND VOIDS RATIO ARE DEFINED VIA SDVINI
C
      VOID=STATEV(1)
      PMAX=STATEV(2)
C
c	write(*,*)'ch1'
C
      CALL KTHETA(STRESS,NTENS,STATEV,NSTATV,PROPS(1),THETA,DETS,
     *AJ,CMN,COS3T,RMTHETA)
C
C     CALCULATE PI (from stress state), C1, C2 AND C3
C
      PRESS=(STRESS(1)+STRESS(2)+STRESS(3))/(-3.0)
      IF(PRESS.LE.0.0)PRESS=1.0
      QVAL=QDS(STRESS,NTENS)
      IF(QVAL.LT.0.001)QVAL=0.001
      ETA=QVAL/PRESS
C
      RNEXP=(RNOV-1.)/RNOV
      C1=RNOV-1.0
      C2=(RNOV)/(1.-RNOV)
      TERM1=((1./(1.-RNOV)-C2*ETA/RMTHETA))
      IF(TERM1.LE.0.0)TERM1=(1./(1.-RNOV)-C2)
      PIPQ=PRESS*(TERM1**RNEXP)
C
C
C EACH ITERATION IS DIVIDED INTO SMALL PORTIONS
C FIRST PORTION IS 1/NDIV, NDIV IS SUGGESTED NO OF DIVISIONS
C
      T=0.0
      DT=1.0/5.0
      ICOUNT=0
C	
110   CONTINUE
      ICOUNT=ICOUNT+1
C
C
C SWITCH BETWEEN CONSTANT AND EXPONENTIAL HARDENING MODULUS
C
	IF(ISWIT.NE.1)HARD=PROPS(5)*EXP(1-ETA/PROPS(1))
	IF(ISWIT.EQ.1)HARD=PROPS(5)
C
C
C CALCULATE SHEAR MODULUS AND BULK MODULUS
C
      GS=PROPS(6)*PRESS**PROPS(7)
      BKS=(2.0*(1+PROPS(8)))/(3.0*(1-2.0*PROPS(8)))*GS
C
C     ELASTIC STIFFNESS MATRIX
C
      DE(1,1)=BKS+4.0/3.0*GS                                            
      DE(1,2)=BKS-2.0/3.0*GS                                                  
      DE(1,3)=DE(1,2)                                                   
      DE(2,1)=DE(1,2)                                                   
      DE(2,2)=DE(1,1)                                                   
      DE(2,3)=DE(1,3)                                                   
      DE(3,1)=DE(1,3)                                                   
      DE(3,2)=DE(2,3)                                                   
      DE(3,3)=DE(1,1)                                                   
      DE(4,4)=GS                                                             
      IF(NDIM.EQ.2)GO TO 111
      DE(5,5)=GS                                                             
      DE(6,6)=GS                                                             
111   CONTINUE
C
C    DECIDE ELASTIC OR PLASTIC (OCS OR NCS)
C    CALCULATE tr(df/ds E de) and if it is ge 0.0 loading lt 0.0 unloading
C    df/ds is A(i,j), E is DE(i,j,k,l) and de is DSTRAN (k,l)
C    HERE CRISP CRITERION IS USED TO DECIDE LOADING OR UNLOADING
C
      CALL KELAPLA(RNOV,C1,C2,PMAX,PIPQ,STATEV,NTENS,STRAN,DSTRAN,
     *PRESS,ETA,RMTHETA, NSTATV,KINC,NPT,NOEL)
C	
C
      IF(STATEV(16).EQ.0.0)THEN
      DO IW=1,NTENS
      DO IX=1,NTENS
      DDSDDE(IW,IX)=DE(IW,IX)
      END DO
      END DO
C
      DO K1=1,NTENS                                                          
      DO K2=1,NTENS	                                                        
      STRESS(K2)=STRESS(K2)+DDSDDE(K2,K1)*DSTRAN(K1)*DT
      END DO                                                                  
      END DO              
C
      DO I5=1,NTENS
      DO I6=1,NTENS
      DSUM(I5,I6)=DSUM(I5,I6)+DDSDDE(I5,I6)
      END DO
      END DO
C
      DVOLS=DSTRAN(1)+DSTRAN(2)+DSTRAN(3)
      DV=DVOLS*DT*(1+VOID)
      VOID=VOID+DV
C
	PRESS=(STRESS(1)+STRESS(2)+STRESS(3))/(-3.0)
      QVAL=QDS(STRESS,NTENS)
      IF(PRESS.LE.0.0)PRESS=1.0
      IF(QVAL.LT.0.001)QVAL=0.001
      ETA=QVAL/PRESS
C
      TERM1=((1./(1.-RNOV)-C2*ETA/RMTHETA))
      IF(TERM1.LE.0.0)TERM1=(1./(1.-RNOV)-C2)
      PIPQ=PRESS*(TERM1**RNEXP)
C
      RAT=(1.0/(1.0-RNOV))**(1.0/C2)
      PMN=PIPQ*RAT
      IF(PMN.GE.PMAX)STATEV(2)=PMN
C
C
C   CHECK WHETEHR INCREMENT HAS BEEN FINISHED
C
      T=T+DT
      IF(T.LT.1)GO TO 110      
C
      GO TO 113
      ENDIF
C
C STRESS STATE IS PLASTIC - PLASTIC CALCULATIONS
C
	CALL KDDSDDE(STRESS,NTENS, STATEV, NSTATV,PRESS,QVAL,PIPQ,
     *C1,C2,RNOV,PROPS,NPROPS,EMIN,EMAX,VOID,RNEXP,DILC,PI_MAX,DE,
     *DDSDDE,HARD,A,PSI,NPT,RMTHETA,ECRIT)
C
	DO I5=1,NTENS
	DO I6=1,NTENS
	DSUM(I5,I6)=DSUM(I5,I6)+DDSDDE(I5,I6)
	END DO
	END DO
C
	IF(T.EQ.0.0)THEN
	DO I7=1,NTENS
	AF(I7)=A(I7)
	END DO
	ENDIF
C        
C UPDATE STRESS STATE, PI, VOID
C
      DO IY=1,NTENS
      STRI1(IY)=0.0
      END DO
C	
      DO K1=1,NTENS                                                          
      DO K2=1,NTENS	                                                       
      STRESS(K2)=STRESS(K2)+DDSDDE(K2,K1)*DSTRAN(K1)*DT
      STRI1(K2)=STRI1(K2)+DDSDDE(K2,K1)*DSTRAN(K1)*DT
      END DO                                                                  
      END DO
C
C
	PRESS=(STRESS(1)+STRESS(2)+STRESS(3))/(-3.0)
	QVAL=QDS(STRESS,NTENS)
	IF(PRESS.LE.0.0)PRESS=1.0
	IF(QVAL.LT.0.001)QVAL=0.001
	ETA=QVAL/PRESS
C
      RNEXP=(RNOV-1.)/RNOV
      C1=RNOV-1.0
      C2=(RNOV)/(1.-RNOV)
	TERM1=((1./(1.-RNOV)-C2*ETA/RMTHETA))
      IF(TERM1.LE.0.0)TERM1=(1./(1.-RNOV)-C2)
      PIPQ=PRESS*(TERM1**RNEXP)
C
C
      DVOLS=DSTRAN(1)+DSTRAN(2)+DSTRAN(3)
      DV=DVOLS*DT*(1+VOID)
      VOID=VOID+DV
C
	T=T+DT
	IF(T.LT.1)GO TO 110
C
113    CONTINUE
C
C
	DO I8=1,NTENS
	DO I9=1,NTENS
	DDSDDE(I8,I9)=DSUM(I8,I9)/ICOUNT
	END DO
	END DO
C
C       DRIFT CORRECTION
C
      PRESS=(STRESS(1)+STRESS(2)+STRESS(3))/(-3.0)
      IF(PRESS.LE.0.0)PRESS=1.0
	IF(PRESS.EQ.0.0)PRESS=1.0
      QVAL=QDS(STRESS,NTENS)
      IF(QVAL.LT.0.001)QVAL=0.001
      ETA=QVAL/PRESS
C
C
      RAT=(1.0/(1.0-RNOV))**(1.0/C2)
      PMN=PIPQ*RAT
      IF(PMN.GE.PMAX)STATEV(2)=PMN
C
C UPDATE STRAN FOR OUTPUT
C
	DO I6=1,NTENS
	STRAN(I6)=STRAN(I6)+DSTRAN(I6)
	END DO
C
C CALCULATE SHEAR STRAIN
C
      EDS1=EDS(STRAN,NTENS)
C
      STATEV(1)=VOID
      STATEV(4)=QVAL
      STATEV(5)=ETA
      STATEV(6)=EDS1
      STATEV(7)=(STRAN(1)+STRAN(2)+STRAN(3))*(-1.0)
      STATEV(8)=PSI
      STATEV(9)=0
      STATEV(10)=RMTHETA
      STATEV(11)=PIPQ
      STATEV(12)=PI_MAX
      STATEV(17)=PRESS
      STATEV(18)=ECRIT
C
      RETURN
      END
C
C CALCULATE DEVIATORIC STRESS
C      
      DOUBLE PRECISION FUNCTION QDS(VARINT,NTENS)                              
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                       
      DIMENSION VARINT(6)                                                       
C                                                                               
      SX=VARINT(1)                                                              
      SY=VARINT(2)                                                              
      SZ=VARINT(3)                                                              
      TXY=VARINT(4)                                                             
      Q2=SX*(SX-SY)+SY*(SY-SZ)+SZ*(SZ-SX)+3.*TXY*TXY                           
      IF(NTENS.EQ.4)GO TO 5
      TYZ=VARINT(5)                                                             
      TZX=VARINT(6) 
      Q2=Q2+3.*TYZ*TYZ+3.*TZX*TZX                                               
5     QDS=SQRT(Q2)                                                              
      RETURN                                                                  
      END
C
C
C CALCULATE DEVIATORIC STRAIN
C
      DOUBLE PRECISION FUNCTION EDS(A,NTENS)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(6)
      EDS2=0.5*((A(1)-A(2))*(A(1)-A(2))+(A(2)-A(3))*(A(2)-A(3))+
     1   (A(3)-A(1))*(A(3)-A(1)))+0.75*A(4)*A(4)
      IF(NTENS.EQ.4)GOTO 10
      EDS2=EDS2+0.75*A(5)*A(5)+0.75*A(6)*A(6)
   10 EDS=2.*SQRT(EDS2)/3.
      RETURN
      END
C
C
	SUBROUTINE KDDSDDE(RLSTRS,NTENS,STATEV,NSTATV,PRESS,QVAL,
     *PI,C1,C2,RNOV,PROPS,NPROPS,EMIN,EMAX,VOID,RNEXP,DILC,PI_MAX,
     *DE,DDSDDE,HARD,A,PSI,NPT,RMTHETA,ECRIT)
C	
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C	
	DIMENSION PROPS(NPROPS),STATEV(NSTATV),DDSDDE(NTENS,NTENS)
	DIMENSION RLSTRS(NTENS)
	DIMENSION SIGMA(3,3),S(3,3),A(6),TOP1(6)
	DIMENSION TOP(6,6),BAS1A(6),DEP(6,6),DE(6,6)
      DIMENSION DJDSG(6),DDSDSG(6),DTDSG(6),DFDSG3(6)
C
	NDIM=3
	IF(NTENS.EQ.4)NDIM=2
C
	TOLEK=1.E-17
	CALL ZERO1(A,6)
	CALL ZERO1(DJDSG,6)
	CALL ZERO1(DDSDSG,6)
	CALL ZERO1(DTDSG,6)
	CALL ZERO1(DFDSG3,6)
	CALL ZERO1(TOP1,6)
	CALL ZERO1(BAS1A,6)
	CALL ZERO2(TOP,6)
	CALL ZERO2(DEP,6)
	CALL ZERO2(SIGMA,3)
	CALL ZERO2(S,3)	
C
	SIGMA(1,1)=-RLSTRS(1)
	SIGMA(2,2)=-RLSTRS(2)
	SIGMA(3,3)=-RLSTRS(3)
	SIGMA(1,2)=-RLSTRS(4)
	SIGMA(2,1)=SIGMA(1,2)
	IF(NDIM.EQ.2)GO TO 122
	SIGMA(1,3)=-RLSTRS(5)
	SIGMA(3,1)=SIGMA(1,3)
	SIGMA(2,3)=-RLSTRS(6)
	SIGMA(3,2)=SIGMA(2,3)
122	CONTINUE	
C
C	CALCULATE Sij
C
	S(1,1)=SIGMA(1,1)-PRESS
	S(2,2)=SIGMA(2,2)-PRESS
	S(3,3)=SIGMA(3,3)-PRESS
	S(1,2)=SIGMA(1,2)
	S(2,1)=S(1,2)
	IF(NDIM.EQ.2)GO TO 123                                          
	S(2,3)=SIGMA(2,3)                                               
	S(3,2)=S(2,3)                                                   
	S(3,1)=SIGMA(3,1)                                               
	S(1,3)=S(3,1)                                                   
123   CONTINUE
C	
C
	CALL KTHETA(RLSTRS,NTENS,STATEV,NSTATV,PROPS(1),THETA,DETS,
     *AJ,CMN,COS3T,RMTHETA)	
C
C			       De a b
C     PLASTIC DMATRIX D = (1 - -------       )De
C			       b De a - c H a
C     LET US CALCULATE a =df/dsigma
C     CALCULATE ADDITIONAL TERMS FOR A VECTOR
C
C***********************************************************************
C
C    CALCULATE DFDP, DFDM, DMDT AND FINALLY DFDSIGMA (A)
C
C***********************************************************************
C
      DFDP=(-RMTHETA/(RNOV*3.))*(1+C1*(1.+C2)*(PRESS/PI)**C2)
      DFDM=-PRESS/RNOV*(1.+C1*(PRESS/PI)**C2)
	DMDT=(-2.0/SQRT(27.0)*CMN*3.0*COS3T*RMTHETA**3/3.0**1.5)/
	1     (2.0/3.0*(CMN-3.0)*RMTHETA+2.0/SQRT(27.0)*CMN*SIN(3.0*THETA)*
	2     RMTHETA**2/SQRT(3.0))
C
C     CALCULATE DJDSG
C
      IF(ABS(AJ).LT.TOLEK)THEN
	DJDSG(1)=1/(2.0)*S(1,1)
	DJDSG(2)=1/(2.0)*S(2,2)
	DJDSG(3)=1/(2.0)*S(3,3)
	DJDSG(4)=1/(2.0)*2.0*S(1,2)
      IF(NDIM.EQ.2)GO TO 125
	DJDSG(5)=1/(2.0)*2.0*S(3,1)
	DJDSG(6)=1/(2.0)*2.0*S(2,3)
	ELSE
	DJDSG(1)=1/(2.0*AJ)*S(1,1)
	DJDSG(2)=1/(2.0*AJ)*S(2,2)
	DJDSG(3)=1/(2.0*AJ)*S(3,3)
	DJDSG(4)=1/(2.0*AJ)*2.0*S(1,2)
	IF(NDIM.EQ.2)GO TO 125
	DJDSG(5)=1/(2.0*AJ)*2.0*S(3,1)
	DJDSG(6)=1/(2.0*AJ)*2.0*S(2,3)
	END IF
C	
C     CALCULATE DDSDSG
C
125	DDSDSG(1)=-1./3.*S(1,1)*S(2,2)-1./3.*S(1,1)*S(3,3)+
	1          2./3.*S(2,2)*S(3,3)-2./3.*S(2,3)**2+
     2          1./3.*S(3,1)**2+1./3.*S(1,2)**2
	DDSDSG(2)=-1./3.*S(1,1)*S(2,2)+2./3.*S(1,1)*S(3,3)-
	1          1./3.*S(2,2)*S(3,3)+1./3.*S(2,3)**2-
     2          2./3.*S(3,1)**2+1./3.*S(1,2)**2
	DDSDSG(3)=2./3.*S(1,1)*S(2,2)-1./3.*S(1,1)*S(3,3)-
	1          1./3.*S(2,2)*S(3,3)+1./3.*S(2,3)**2+
     2          1./3.*S(3,1)**2-2./3.*S(1,2)**2
	DDSDSG(4)=-2.0*S(3,3)*S(1,2)+2.0*S(2,3)*S(3,1)
	IF(NDIM.EQ.2)GO TO 126
	DDSDSG(5)=-2.0*S(2,2)*S(3,1)+2.0*S(1,2)*S(2,3)
	DDSDSG(6)=-2.0*S(1,1)*S(2,3)+2.0*S(1,2)*S(3,1)
C
C	CALCULATE DTDSG
C
C
126	IF(ABS(AJ).LT.TOLEK.OR.ABS(COS3T).LT.TOLEK)THEN
	DO INN=1,NTENS
	DTDSG(INN)=0.0
	END DO
	ELSE
	DTDSG(1)=SQRT(3.0)/(2.0*COS3T*AJ**3)*
	1         (3.0*DETS/AJ*DJDSG(1)-DDSDSG(1))
	DTDSG(2)=SQRT(3.0)/(2.0*COS3T*AJ**3)*
	1         (3.0*DETS/AJ*DJDSG(2)-DDSDSG(2))
	DTDSG(3)=SQRT(3.0)/(2.0*COS3T*AJ**3)*
	1         (3.0*DETS/AJ*DJDSG(3)-DDSDSG(3))
	DTDSG(4)=SQRT(3.0)/(2.0*COS3T*AJ**3)*
	1         (3.0*DETS/AJ*DJDSG(4)-DDSDSG(4))
      IF(NDIM.EQ.2)GO TO 127
	DTDSG(5)=SQRT(3.0)/(2.0*COS3T*AJ**3)*
	1         (3.0*DETS/AJ*DJDSG(5)-DDSDSG(5))
	DTDSG(6)=SQRT(3.0)/(2.0*COS3T*AJ**3)*
	1         (3.0*DETS/AJ*DJDSG(6)-DDSDSG(6))
	END IF	
C
C	CALCULATE DFDSG3
C
127	DFDSG3(1)=DFDM*DMDT*DTDSG(1)
	DFDSG3(2)=DFDM*DMDT*DTDSG(2)
	DFDSG3(3)=DFDM*DMDT*DTDSG(3)
	DFDSG3(4)=DFDM*DMDT*DTDSG(4)
	IF(NDIM.EQ.2)GO TO 128
	DFDSG3(5)=DFDM*DMDT*DTDSG(5)
	DFDSG3(6)=DFDM*DMDT*DTDSG(6)
C
C	CALCULATE DFDSG=A
C
128   A(1)=DFDP+(1./(2.*QVAL))*(-3.*RLSTRS(1)-3.*PRESS)+DFDSG3(1)
C      
      A(2)=DFDP+(1./(2.*QVAL))*(-3.*RLSTRS(2)-3.*PRESS)+DFDSG3(2)
C
      A(3)=DFDP+(1./(2.*QVAL))*(-3.*RLSTRS(3)-3.*PRESS)+DFDSG3(3)
C
      A(4)=(1./QVAL)*(-3*RLSTRS(4))+DFDSG3(4)
C
      IF(NDIM.EQ.2)GO TO 129	
      A(5)=(1./QVAL)*(-3*RLSTRS(5))+DFDSG3(5)
C
      A(6)=(1./QVAL)*(-3*RLSTRS(6))+DFDSG3(6)
C
C***********************************************************************
C
C    CALCULATE THE MATRIX
C
C***********************************************************************
C
C     EVALUATE De a 
C
129   DO IK=1,NTENS
      TOP1(IK)=0.0
      DO IL=1,NTENS
      TOP1(IK)=TOP1(IK)+DE(IK,IL)*A(IL)
      END DO
      END DO
C
C
C     Evaluate (De a) bt
C
      DO IM=1,NTENS
      DO IN=1,NTENS
      TOP(IM,IN)=TOP1(IM)*A(IN)
      END DO
      END DO
C
C     Evaluate b De
C
      DO IO=1,NTENS
      BAS1A(IO)=0.0
      DO IP=1,NTENS
      BAS1A(IO)=BAS1A(IO)+A(IP)*DE(IP,IO)
      END DO
      END DO
C
C     Evaluate (b De) a
C
      BAS1=0.0
      DO IQ=1,NTENS
      BAS1=BAS1+BAS1A(IQ)*A(IQ)
      END DO
C
C    NOW FORMULATE C AND H (C)= (RC1,0)
C      
C    FIRST CALCULATE PI_MAX AND DILATION
C
C    THIS IS WHERE THE EMAX AND EMIN INTRODUCED
C
      ECRIT=EMAX-(EMAX-EMIN)/(9.9035-LOG(PRESS))
      PSI=VOID-ECRIT
C
      ECRIT=EMAX-(EMAX-EMIN)/(9.9035-LOG(PRESS))
      PSI=VOID-ECRIT
      ETA=QVAL/PRESS
      IF(RMTHETA.LT.TOLEK.OR.(1./(1.-RNOV)-C2*ETA/RMTHETA).LT.TOLEK)THEN
	PI_BY_P=(1./(1.-RNOV)-C2)**(RNEXP)
	ELSE
	PI_BY_P=(1./(1.-RNOV)-C2*ETA/RMTHETA)**(RNEXP)
	END IF
	IF(PI_BY_P.LE.0.0)PI_BY_P=(1./(1-RNOV)-C2)**(RNEXP)
	PI_BY=PI_BY_P*PRESS
	EI=EMAX-(EMAX-EMIN)/(9.9035-LOG(PI_BY))
	PSI_I=VOID-EI
      DMAX=DILC*PSI_I*RMTHETA/PROPS(1)
      PI_MAX=(abs((1.+DMAX*RNOV/PROPS(1))))**(RNEXP)*PRESS
      DIL=(RMTHETA-ETA)/(1.-RNOV)
C
C      
      RC1=C1*C2*(RMTHETA/RNOV)*(PRESS/PI)**(1+C2)*HARD*(PI_MAX-PI)*
	1     RMTHETA/PROPS(1)
C
      BAS2=RC1*SQRT(2./3.)*SQRT((A(1)-(A(1)+A(2)+A(3))/3)**2+
	1     (A(2)-(A(1)+A(2)+A(3))/3)**2+(A(3)-(A(1)+A(2)+A(3))/3)**2+
     2     2*A(4)**2+2*A(5)**2+2*A(6)**2)   
C 
C     NOW FORMULATE D-MATRIX
C
      DO IW=1,NTENS
      DO IX=1,NTENS
      TOP(IW,IX)=TOP(IW,IX)/(BAS1-BAS2)
      END DO
      END DO
C
C     MULTIPLY DE WITH TOP AND CALL IT DEP
C
      DO IY=1,NTENS
      DO IZ=1,NTENS
      DEP(IY,IZ)=0.0      
      DO I1=1,NTENS 
      DEP(IY,IZ)=DEP(IY,IZ)+TOP(IY,I1)*DE(I1,IZ)
      END DO
      END DO
      END DO
C                 
      DO IW=1,NTENS
      DO IX=1,NTENS
      DDSDDE(IW,IX)=DE(IW,IX)-DEP(IW,IX)
      END DO
      END DO
C
      RETURN
      END	

C
C***********************************************************************
C
	SUBROUTINE KPLAS(DSTRAN,DE,NTENS,PLAS,STRI1,DT)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                                      
	DIMENSION DSTRAN(NTENS),STRI1(6)
C
	DIMENSION A(6,6),ELAS(6),PLAS(6),DE(6,6)
C
	CALL ZERO1(ELAS,NTENS)
	DO IU=1,NTENS
	DO IY=1,NTENS
	A(IU,IY)=DE(IU,IY)
	END DO
	END DO
C
	NMAX=6
	CALL KINV(A,NTENS,NMAX)
	DO IU=1,NTENS
	DO IY=1,NTENS
	ELAS(IU)=ELAS(IU)+A(IU,IY)*STRI1(IY)
	END DO
	END DO
C	
	DO IZ=1,NTENS
	PLAS(IZ)=DSTRAN(IZ)*DT-ELAS(IZ)
	END DO		
C	
	RETURN
	END 		
C
C***********************************************************************
C
C	INVERSE OF A MATRIX
C
C***********************************************************************
C
	SUBROUTINE KINV(A,N,NMAX)
C
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
	DIMENSION A(NMAX,NMAX),IS(60),JS(60)
C
        DO 100 K=1,N
          D=0.
          DO 10 I=K,N
            DO 10 J=K,N
              IF(ABS(A(I,J)).GT.D)THEN
                D=ABS(A(I,J))
                IS(K)=I
                JS(K)=J
              ENDIF
10        CONTINUE
          IF(D+1.0.EQ.1.0)THEN
            WRITE(6,20)
            DO IY=1,N
            WRITE(6,21)(A(IY,IZ),IZ=1,N)
            END DO
            STOP
          ENDIF
20        FORMAT(1X,'ERR * * NOT INV')
21	  FORMAT(4(1X,F12.5))
          DO 30 J=1,N
            T=A(K,J)
            A(K,J)=A(IS(K),J)
            A(IS(K),J)=T
30        CONTINUE
          DO 40 I=1,N
            T=A(I,K)
            A(I,K)=A(I,JS(K))
            A(I,JS(K))=T
40        CONTINUE
          A(K,K)=1./A(K,K)
          DO 50 J=1,N
            IF(J.NE.K)THEN
              A(K,J)=A(K,J)*A(K,K)
            ENDIF
50        CONTINUE
          DO 70 I=1,N
            IF(I.NE.K)THEN
              DO 60 J=1,N
                IF(J.NE.K)THEN
                  A(I,J)=A(I,J)-A(I,K)*A(K,J)
                ENDIF
60            CONTINUE
            ENDIF
70        CONTINUE
          DO 80 I=1,N
            IF(I.NE.K)THEN
              A(I,K)=-A(I,K)*A(K,K)
            ENDIF
80        CONTINUE
100     CONTINUE
C
        DO 130 K=N,1,-1
          DO 110 J=1,N
            T=A(K,J)
            A(K,J)=A(JS(K),J)
            A(JS(K),J)=T
110       CONTINUE
          DO 120 I=1,N
            T=A(I,K)
            A(I,K)=A(I,IS(K))
            A(I,IS(K))=T
120       CONTINUE
130     CONTINUE
        RETURN
        END
C
C***********************************************************************
C
	SUBROUTINE ZERO2(A,N)
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C	
	DIMENSION A(N,N)
	DO I=1,N
	DO J=1,N
	A(I,J)=0.0
	END DO
	END DO
	RETURN
	END
C
	SUBROUTINE ZERO1(A,N)
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C	
	DIMENSION A(N)
	DO I=1,N
	A(I)=0.0
	END DO
	RETURN
	END
C
	SUBROUTINE KNEWDT(RK1,TOL,DT)
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C	
	DTOLD=DT
	IF(RK1.LE.1E-5)RK1=1E-5
	Q=0.8*SQRT((TOL/RK1))
	DT=Q*DTOLD
	IF(DT.LE.0.1*DTOLD)DT=0.1*DTOLD
	IF(DT.GE.2.5*DTOLD)DT=2.5*DTOLD
	RETURN
	END
C
C
C***********************************************************************
C
C    DECIDE WHETHER STRESS STATE IS INSIDE OR ON THE YIELD SURFACE
C
C***********************************************************************
C
	SUBROUTINE KELAPLA(RNOV,C1,C2,PMAX,PIPQ,STATEV,NTENS,STRAN,
     *DSTRAN,PRESS,ETA,RMTHETA,NSTATV,KINC,NPT,NOEL)
C
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C	
	DIMENSION STATEV(NSTATV),STRAN(NTENS),DSTRAN(NTENS)
	DIMENSION STRAN2(6)
C
	CALL ZERO1(STRAN2,NTENS)
C
      RAT=(1.0/(1.0-RNOV))**(1.0/C2)
      PIM=PMAX/RAT
      PMN=PIPQ*RAT
C            		
      ED1=STATEV(6)
      DO IY=1,NTENS
      STRAN2(IY)=STRAN(IY)+DSTRAN(IY)
      END DO
      ED2=EDS(STRAN2,NTENS)
C
      IF(PRESS.GE.PIPQ)THEN
	      IF(PMN.LT.0.99*PMAX)THEN
	      STATEV(16)=0.0
	      STATEV(3)=0.0
	      ENDIF
         IF(PMN.GE.0.99*PMAX)THEN
         STATEV(16)=1.0
         STATEV(3)=1.0
         ENDIF
      ENDIF
C      
      IF(PRESS.LT.PIPQ.AND.ETA.LE.RMTHETA)THEN
	      IF(PMN.LT.0.99*PMAX)THEN
	      STATEV(16)=0.0
	      STATEV(3)=2.0
	      ENDIF
         IF(PMN.GE.0.99*PMAX)THEN
         STATEV(16)=1.0
         STATEV(3)=4.0
         ENDIF
      ENDIF
C
      IF(PRESS.LT.PIPQ.AND.ETA.GT.RMTHETA)THEN
	      IF(PMN.GE.PMAX)THEN
     	      STATEV(16)=1.0
	      STATEV(3)=4.0
	      ENDIF
		    IF(PMN.LT.PMAX)THEN
	      		IF(ED2.GE.ED1)THEN
	    	    STATEV(16)=1.0
 	     	    STATEV(3)=4.0
		        ENDIF
	                  IF(ED2.LT.ED1)THEN
      	              STATEV(3)=3.0	
      	              STATEV(16)=0.0
      	              ENDIF
		    ENDIF
      ENDIF      
	RETURN
	END
C
C***********************************************************************
C
C    CALCULATE THETA AND RMTHETA
C
C***********************************************************************
C
	SUBROUTINE KTHETA(RLSTRS,NTENS,STATEV,NSTATV,PROPS1,THETA,DETS,
     *AJ,CMN,COS3T,RMTHETA)
C
	IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C	
	DIMENSION STATEV(NSTATV),RLSTRS(NTENS)
	DIMENSION SIGMA(3,3),S(3,3)
C
	CALL ZERO2(SIGMA,3)
	CALL ZERO2(S,3)	
C
	NDIM=3
	IF(NTENS.EQ.4)NDIM=2
      PRESS=(RLSTRS(1)+RLSTRS(2)+RLSTRS(3))/(-3.0)
C
	SIGMA(1,1)=-RLSTRS(1)
	SIGMA(2,2)=-RLSTRS(2)
	SIGMA(3,3)=-RLSTRS(3)
	SIGMA(1,2)=-RLSTRS(4)
	SIGMA(2,1)=SIGMA(1,2)
	IF(NDIM.EQ.2)GO TO 122
	SIGMA(1,3)=-RLSTRS(5)
	SIGMA(3,1)=SIGMA(1,3)
	SIGMA(2,3)=-RLSTRS(6)
	SIGMA(3,2)=SIGMA(2,3)
122	CONTINUE	
C
C	CALCULATE Sij
C
	S(1,1)=SIGMA(1,1)-PRESS
	S(2,2)=SIGMA(2,2)-PRESS
	S(3,3)=SIGMA(3,3)-PRESS
	S(1,2)=SIGMA(1,2)
	S(2,1)=S(1,2)
	IF(NDIM.EQ.2)GO TO 123                                          
	S(2,3)=SIGMA(2,3)                                               
	S(3,2)=S(2,3)                                                   
	S(3,1)=SIGMA(3,1)                                               
	S(1,3)=S(3,1)                                                   
123     CONTINUE
C
C***********************************************************************
C
C    MATSUOKA AND NAKAI'S FAILURE CRITERIA
C    CALCULATING LODE ANGLE
C
C***********************************************************************
C
C	CALCULATE THE VALUE OF DETS
C
	DETS=S(1,1)*S(2,2)*S(3,3)-S(1,1)*S(2,3)**2-S(2,2)*S(3,1)**2-
	1     S(3,3)*S(1,2)**2+2*S(1,2)*S(2,3)*S(3,1)
C
C	CALCULATE THE VALUE OF J AS AJ
C
	AJ=1.0/SQRT(6.0)*SQRT((SIGMA(1,1)-SIGMA(2,2))**2+(SIGMA(2,2)-
	1   SIGMA(3,3))**2+(SIGMA(3,3)-SIGMA(1,1))**2+SIGMA(1,2)**2+
     2   SIGMA(2,3)**2+SIGMA(3,1)**2)
	IF(AJ.GT.-1.E-8.AND.AJ.LT.1.E-8)AJ=0.0
C
C	CALCULATE DETS/AJ**3
C
	IF(AJ.EQ.0.0)THEN
	DETSAJ3=0.0
	ELSE
	DETSAJ3=DETS/AJ**3
	END IF
C
C	CALCULATE THETAV
C
	THETAV=3.0*SQRT(3.0)/2.0*DETSAJ3
	IF(THETAV.GT.0.99)THETAV=1.0
	IF(THETAV.LT.-0.99)THETAV=-1.0
C
C	CALCULATE THE VALUE OF THETA
C
	THETA=-1.0/3.0*ASIN(THETAV)  
	IF(THETA.GT.0.523598)THETA=0.523598
	IF(THETA.LT.-0.523598)THETA=-0.523598 
      STATEV(13)=THETA
C
C	CALCULATE COS3THETA
C
	COS3T=COS(3.0*THETA)
	IF(COS3T.GT.-1.E-8.AND.COS3T.LT.1.E-8)COS3T=0.0
C
C	CALCULATE MJ
C
	RMJ=PROPS1/SQRT(3.0)
C
C	CALCULATE CMN
C
	CMN=(9.0-3.0*RMJ**2)/(2.0*SQRT(3.0)/9.0*RMJ**3-
	1    RMJ**2+1.0)
C
C	CALCULATE AJ2ETHA USING NEWTON'S METHOD
C
	AJ2ETHA=0.1
	TOLE=1.E-5
	NN=0
195	IF(NN.LT.50)THEN
	FX=(CMN-3.0)*AJ2ETHA+2.0/SQRT(27.0)*CMN*SIN(3.0*THETA)*
	1   AJ2ETHA**1.5-(CMN-9.0)
	FDX=(CMN-3.0)+2.0/SQRT(27.0)*CMN*3.0/2.0*SIN(3.0*THETA)*
	1    AJ2ETHA**0.5
	AJ2ETHA1=AJ2ETHA-(FX/FDX)
	END IF
	IF(ABS(AJ2ETHA1-AJ2ETHA).GT.TOLE.AND.NN.LT.50)THEN
	AJ2ETHA=AJ2ETHA1
	NN=NN+1
	GO TO 195
	ELSE
	AJ2ETHA=AJ2ETHA1
	END IF
C
C	OBTAIN THE VALUE OF RMTHETA
C
	RMTHETA=SQRT(3.0*AJ2ETHA)
	STATEV(14)=RMTHETA
	RETURN
	END
C	
C	
C***********************************************************************
C
C    SUBROUTINE TO INPUT INITIAL VOID RATIO
C
C***********************************************************************
C
      SUBROUTINE SDVINI(STATEV,COORDS,NSTATV,NCRDS,NOEL,NPT,LAYER,KSPT)         
C
      INCLUDE 'ABA_PARAM.INC'                                                   
      DOUBLE PRECISION STATEV(NSTATV),COORDS(NCRDS)                             
C                                                                               
C     Y1 IS TOP OF MESH AND Y2 IS BASE OF THE MESH                              
C
C     E1 AND E2 ARE VOID RATIOS  (STATEV 1)
C     PMAX1 PMAX2 ARE PMAXs at TOP AND BASE (STATEV 2)
C                                                                               
      STATEV(1)=0.55
C      E1=0.669                                                               
C      E2=0.669                                                                 
C
C      PMAX1=300.0                                                               
C      PMAX2=300.0                                                               
C      Y1=1.26                                                                 
C      Y2=-1.0	                                                                
C                                                                               
C      YC=COORDS(2)                                                              
C      STATEV(1)=((E2-E1)/(Y1-Y2))*((Y1-Y2)-YC)+E1                               
C      STATEV(2)=((PMAX2-PMAX1)/(Y1-Y2))*((Y1-Y2)-YC)+PMAX1                      
C                                                                               
      RETURN                                                                    
      END   
C
C***********************************************************************
C                                                                    