       module kvisual
       implicit none
       double precision sigout(32,3,4)
       save
       end module    

	   SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     1 PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
     2 KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,NPREDF,
     3 LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,PERIOD)
C
	   use kvisual
       INCLUDE 'ABA_PARAM.INC'
C
	  
       dimension RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),PROPS(*),
     1 SVARS(NSVARS),ENERGY(8),COORDS(MCRD,NNODE),U(NDOFEL),
     2 DU(MLVARX,*),V(NDOFEL),A(NDOFEL),TIME(2),PARAMS(*),
     3 JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),DDLMAG(MDLOAD,*),
     4 PREDEF(2,NPREDF,NNODE),LFLAGS(*),JPROPS(*)

	   PARAMETER (maxDOF=2 , ndof=2 , ndim=2 ,ninpt=3,nelem=32,nsdv=4,nsvint=4,ntens=4)
	   PARAMETER (gaussCoordWieght = 0.3333333333333,ro=7850e-12)
       double precision dshape(ndim,3),xjac(maxDOF,maxDOF),xjaci(maxDOF,maxDOF),
     1 bmat(3,2*nnode),bmatT(2*nnode,3),C(3,3),iindex(3),kmat(NDOFEL,NDOFEL),
     2 deps(3),stress(3),SIG(4),SDV(4),statevLocal(nsvint),shapef(2,6),
     4 db(3,2*nnode),XYDerivation(ndim,3),force(6),dforce(6),kdamp(NDOFEL,NDOFEL),
     5 SVARS_local(ninpt*nsvint),kmass(NDOFEL,NDOFEL)
	 
	 
		rhs(:,1)=0.d0
		amatrx=0.d0
		C=0.d0				
		force=0.d0
		kmat=0.0
		
		!Elasticity matrix
		E = props(1)
		xnu = props(2)
		xa = E*(1.-xnu)/(1.d0 + xnu)/(1.0-2*xnu)
		xb = xnu/(1.0-xnu)
		xc = (1.d0-2.0*xnu)/2.d0/(1.0-xnu)
		
	    C(1, 1) = xa
	    C(1, 2) = xa*xb
	    C(2, 1) = C(1, 2)
	    C(2, 2) = C(1, 1)
	    C(3, 3) = xa*xc

        do kintk = 1, ninpt	 
			deps=0.	
			statevLocal=0.d0
			call KshapeFunc(kintk,ninpt,nnode,ndim,shapef,dshape)
			
			djac = 1.d0
			call Kjacobian(ndim,COORDS,nnode,djac,dshape,xjaci)
			
			call KBmatrix(xjaci,dshape,nnode,ndim,bmat,bmatT)

            SVARS_local(1:ninpt*nsvint)=SVARS(1:ninpt*nsvint)
			call Kstatevar(kintk,nsvint,SVARS_local,ninpt*nsvint,statevLocal,1)

			SIG=0.d0
			stress=0.d0
			do k1=1,nsdv
				SIG(k1) = statevLocal(k1)
			end do
			
			stress(1)=SIG(1)
			stress(2)=SIG(2)
			stress(3)=SIG(4)
			   
			call Kstraininc(ntens,ndof,ndim,nnode,MLVARX,bmat,DU,deps)		

			stress=stress+matmul(C,deps)
			force=force+matmul(bmatT,stress*djac*gaussCoordWieght)
			
			db=0.0
			db=matmul(C,bmat)
			kmat=kmat+matmul(bmatT,db*gaussCoordWieght*djac)
			
			kmass=kmass+matmul(Transpose(shapef),shapef*gaussCoordWieght*djac)
			
			SIG(1)=stress(1)
			SIG(2)=stress(2)
			SIG(3)=xnu*(stress(1)+stress(2))
			SIG(4)=stress(3)
			
			do k1=1,nsdv
				statevLocal(k1) = SIG(k1)
			end do

			call Kstatevar(kintk,nsvint,SVARS_local,ninpt*nsvint,statevLocal,0)
            SVARS(1:ninpt*nsvint)=SVARS_local(1:ninpt*nsvint)
			
			do k1=1,nsdv
				sigout(jelem,kintk,k1)=SIG(k1)
			end do
		
		end do 
       
        kmat=0.5*kmat
		kmass=ro*0.5*kmass		

		IF (LFLAGS(3).EQ.1) THEN
            IF (LFLAGS(1).EQ.11 .OR. LFLAGS(1).EQ.12) THEN
			ALPHA = PARAMS(1)
			BETA  = PARAMS(2)
			GAMMA = PARAMS(3)
			
			DADU = 1.0/(BETA*DTIME**2)
			DVDU = GAMMA/(BETA*DTIME)
			
			massPropDamp=0.0006
			stiffPropDamp=0.00026

			kdamp=massPropDamp*kmass + stiffPropDamp*kmat

			AMATRX=DADU*kmass + DVDU*(1.0+ALPHA)*kdamp + (1.0+ALPHA)*kmat
			
			RHS(1:6, 1) = -matmul(kmass,A)
     1                  -(1.0+ALPHA)*(matmul(kdamp,V)+matmul(kmat,U))
     2                  + ALPHA*SVARS(ninpt*nsvint+1:NSVARS)			
			SVARS(ninpt*nsvint+1:NSVARS)=matmul(kdamp,V)+matmul(kmat,U)
            END IF
		END IF		

      RETURN
      END
	  
	  
	  
******************************************************
** shape functions                                   *
******************************************************	
		subroutine KshapeFunc(kintk,ninpt,nnode,ndim,shapef,dshape)
		
		include 'aba_param.inc'
		
		parameter (gaussCoord1 = 0.1666666666667,gaussCoord2 = 0.6666666666667)
		double precision shapef(2,6),dshape(ndim,3),coord23(2,3)
		
		data coord23  /gaussCoord1 , gaussCoord1 ,
	2				   gaussCoord2 , gaussCoord1 ,
	3				   gaussCoord1 , gaussCoord2 /
			
		shapef=0.
		dshape=0.
		r=coord23(1,kintk)
		s=coord23(2,kintk)
		
		shapef(1,1)=1.0-r-s
		shapef(1,3)=r
		shapef(1,5)=s
		
		shapef(2,2)=1.0-r-s
		shapef(2,4)=r
		shapef(2,6)=s
	
		! derivation to r
		
		dshape(1,1) = -1.
		dshape(1,2) = 1.
        dshape(1,3) = 0.
	
	
	    ! derivation to s
		
		dshape(2,1) = -1.
		dshape(2,2) = 0.
        dshape(2,3) = 1.
	
	
		return
		end
		
******************************************************
** Jacobian                                          *
******************************************************
		subroutine Kjacobian(ndim,coords,nnode,djac,dshape,xjaci)
		
		include 'aba_param.inc'
		
		dimension xjac(ndim,ndim),xjaci(ndim,ndim)
		dimension coords(2,nnode),dshape(ndim,nnode)

		xjac=0.d0
		xjaci=0.d0

        xjac=matmul(dshape,transpose(coords))
		
		djac=xjac(1,1)*xjac(2,2)-xjac(1,2)*xjac(2,1)
		 
		if (djac .gt. 0.) then
			xjaci(1,1)=xjac(2,2)/djac
			xjaci(2,2)=xjac(1,1)/djac
			xjaci(1,2)=-xjac(1,2)/djac
			xjaci(2,1)=-xjac(2,1)/djac			
		end if
		
		return
		end
******************************************************
** B-matrix                                          *
******************************************************
		subroutine KBmatrix(xjaci,dshape,nnode,ndim,bmat,bmatT)
		
		include 'aba_param.inc'
		
		dimension xjaci(ndim,ndim),dshape(ndim,3),XYDerivation(ndim,3)
		dimension bmat(3,2*nnode),bmatT(2*nnode,3)

		XYDerivation=0.d0

		bmat=0.d0
		bmatT=0.d0
		
		XYDerivation=matmul(xjaci,dshape)
		bmat(1,1)=XYDerivation(1,1)
		bmat(1,3)=XYDerivation(1,2)
		bmat(1,5)=XYDerivation(1,3)
		
		bmat(2,2)=XYDerivation(2,1)
		bmat(2,4)=XYDerivation(2,2)
		bmat(2,6)=XYDerivation(2,3)
		
		bmat(3,1)=XYDerivation(2,1)
		bmat(3,2)=XYDerivation(1,1)
		bmat(3,3)=XYDerivation(2,2)
		bmat(3,4)=XYDerivation(1,2)
		bmat(3,5)=XYDerivation(2,3)
		bmat(3,6)=XYDerivation(1,3)
		
        bmatT=transpose(bmat)
		
		return
		end
		
c*****************************************************************
		  subroutine Kstraininc(ntens,ndof,ndim,nnode,mlvarx,bmat,du,deps)

		  include 'aba_param.inc'
		  dimension deps(3),bmat(3,2*nnode),du(mlvarx,1),xdu(2) 
		  
			
		  deps = 0.d0
		  xdu=0.d0
			
		  do nodi = 1, nnode
			   
		   incr_row = (nodi - 1)*ndof
			
		   do i = 1, ndof         
				xdu(i)= du(i + incr_row,1)
		   end do
			dNidx = bmat(1,1 + incr_row)
			dNidy = bmat(2,1 + incr_row + 1)
			
		   deps(1) = deps(1) + dNidx*xdu(1)
		   deps(2) = deps(2) + dNidy*xdu(2)
		   deps(3) = deps(3) + dNidy*xdu(1) + dNidx*xdu(2)
					   
		  end do

		  return
		  end


		
************************************************************************		
		 subroutine Kstatevar(npt,nsvint,SVARS_local,n,statev_ip,icopy)
	

		  include 'aba_param.inc'

		  dimension SVARS_local(n),statev_ip(nsvint)

		  isvinc= (npt-1)*nsvint     

		  if (icopy .eq. 1) then

			do i = 1, nsvint
			  statev_ip(i)=SVARS_local(i+isvinc)
			end do
c
		  else
c
			do i = 1, nsvint
			  SVARS_local(i+isvinc)=statev_ip(i)
			end do
		  end if

		  return
		  end

************************************************************************
      subroutine umat(stress,statev,ddsdde,sse,spd,scd,
     1 rpl,ddsddt,drplde,drpldt,
     2 stran,dstran,time,dtime,temp,dtemp,predef,dpred,cmname,
     3 ndi,nshr,ntens,nstatv,props,nprops,coords,drot,pnewdt,
     4 celent,dfgrd0,dfgrd1,noel,npt,layer,kspt,kstep,kinc)
C
		use kvisual
      include 'aba_param.inc'
C
      character*80 cmname
      dimension stress(ntens),statev(nstatv),
     1 ddsdde(ntens,ntens),
     2 ddsddt(ntens),drplde(ntens),
     3 stran(ntens),dstran(ntens),time(2),predef(1),dpred(1),
     4 props(nprops),coords(3),drot(3,3),dfgrd0(3,3),dfgrd1(3,3)
	  PARAMETER (maxDOF=2 , ndof=2 , ndim=2 ,ninpt=3,nelem=32,nsdv=4)
	  integer kelem
      
	   ddsdde=0.d0

      kelem=int(noel-nelem)
	 
      do k1 = 1, nsdv
        statev(k1) = sigout(kelem,npt,k1)
      end do 

      return
      end 