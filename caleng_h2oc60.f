      program test
      implicit double precision(a-h,o-z)
      dimension com(3), Eulang(3) 
      Pi=dacos(-1.0d0)
      
      do kk=1,6
      com(1)=0.0+(kk-1)*0.1d0
      do jj=1,6
      com(2)=0.0+(jj-1)*0.1d0
      do ii=1,6
      com(3)=0.0+(ii-1)*0.1d0


      do k=1,7
      Eulang(1)=0.0d0+(k-1)*60.d0*Pi/180.d0
      do j=1,4
      Eulang(2)=0.0d0+(j-1)*60.d0*Pi/180.d0
      do i=1,7
      Eulang(3)=0.0d0+(i-1)*60.d0*Pi/180.d0
 
      call caleng(com, energy, Eulang)!energy is in K
      
      write(33,70)com,Eulang(1:3)*180.d0/Pi, energy
      
      enddo
      enddo
      enddo
      enddo
      enddo
      enddo
  70  format(3F6.2,'  phi =',F5.1,'  the =',F5.1,'  chi =',F5.1,
     + '  energy =', F18.4)
      end ! end program test


c     ==================================================
      subroutine caleng(com, E_H2OC60, Eulang)
c     ==================================================
c     Initially:
c     __________________________________________________
c     this subroutine calculates LJ potential between
c     rigid water and C60 given the coordinates
c     of centre of mass of water within the cage and 
c     the Euler angles describing the rotation of H2O
c     input: com, rotmat, 
c        ROwf, RH1wf, RH2wf, RCwf(C60)
c     output: EH2OC60
c      phi=Eulang(1)
c      theta=Eulang(2)
c      chi=Eulang(3)
c     __________________________________________________
c     rotmat, computed within the code with Eulang
c     ROwf, RH1wf, RH2wf, RCwf put as data   
      implicit double precision(a-h,o-z)
      dimension ROwf(3), RH1wf(3), RH2wf(3), 
     +          com(3), Eulang(3),rotmat(3,3),
     +          RO_sf(3),RH1_sf(3), RH2_sf(3),
     +          RC_sf(60,3), RCwf(60,3)
c  L-J parameters and conversion factor for energy. 
c  A.B. Farimani, Y. Wu, and N.R. Aluru, 
c  Phys. Chem. Chem. Phys. 15, 17993 (2013).
      parameter(epsco=0.1039d0,sigco=3.372d0,
     +          epsch=0.0256d0,sigch=2.640d0,
     +          kcal2k=503.218978939d0)
c H2O coordinates 
      data ROwf/0.d0,0.d0,-0.06563807d0/,
     +     RH1wf/0.7575d0,0.d0,0.52086193d0/,
     +     RH2wf/-0.7575d0,0.d0,0.52086193d0/
c C60 coordinates
      data RCwf(1,:)/1.00131784714278d0,0.7275d0,-3.325555308452120d0/,
     + RCwf(2,:)/1.958591472830480d0,1.4230d0,-2.594263331494400d0/,   
     + RCwf(3,:)/2.959909319973260d0,0.6955d0,-1.829324563361060d0/,   
     + RCwf(4,:)/2.959909319973260d0,-0.6955d0,-1.829324563361060d0/,  
     + RCwf(5,:)/1.958591472830480d0,-1.4230d0,-2.594263331494400d0/,  
     + RCwf(6,:)/1.001317847142780d0,-0.7275d0,-3.325555308452120d0/,  
     + RCwf(7,:)/-0.382469384066670d0,-1.17712d0,-3.32555530845212d0/, 
     + RCwf(8,:)/-1.237696926152220d0,0.00000d0,-3.325555308452120d0/, 
     + RCwf(9,:)/-0.382469384066670d0,1.17712d0,-3.325555308452120d0/, 
     + RCwf(10,:)/-0.748115372545527d0,2.30246d0,-2.59426333149440d0/, 
     + RCwf(11,:)/0.253202474597252d0,3.02996d0,-1.829324563361060d0/, 
     + RCwf(12,:)/1.576122088763810d0,2.60012d0,-1.829324563361060d0/, 
     + RCwf(13,:)/2.341060856897150d0,2.60012d0,-0.591627637208844d0/, 
     + RCwf(14,:)/3.196288398982700d0,1.4230d0,-0.591627637208844d0/,  
     + RCwf(15,:)/3.422270047712680d0,0.7275d0,0.591627637208844d0/,   
     + RCwf(16,:)/3.422270047712680d0,-0.7275d0,0.591627637208844d0/,  
     + RCwf(17,:)/3.196288398982700d0,-1.4230d0,-0.591627637208844d0/, 
     + RCwf(18,:)/2.341060856897150d0,-2.60012d0,-0.591627637208844d0/,
     + RCwf(19,:)/1.576122088763810d0,-2.60012d0,-1.82932456336106d0/, 
     + RCwf(20,:)/0.253202474597251d0,-3.02996d0,-1.82932456336106d0/, 
     + RCwf(21,:)/-0.748115372545527d0,-2.30246d0,-2.5942633314944d0/, 
     + RCwf(22,:)/-1.985812298697750d0,-2.30246d0,-1.8293245633616d0/, 
     + RCwf(23,:)/-2.803421584636580d0,-1.17712d0,-1.8293245633616d0/, 
     + RCwf(24,:)/-2.420952200569910d0,0.00000d0,-2.59426333149440d0/, 
     + RCwf(25,:)/-2.803421584636580d0,1.17712d0,-1.829324563361060d0/,
     + RCwf(26,:)/-1.985812298697740d0,2.30246d0,-1.829324563361060d0/,
     + RCwf(27,:)/-1.749433219688310d0,3.02996d0,-0.591627637208844d0/,
     + RCwf(28,:)/-0.365645988478857d0,3.47958d0,-0.591627637208844d0/,
     + RCwf(29,:)/0.365645988478858d0,3.47958d0,0.591627637208844d0/,  
     + RCwf(30,:)/1.749433219688310d0,3.02996d0,0.591627637208844d0/,  
     + RCwf(31,:)/1.985812298697750d0,2.30246d0,1.829324563361060d0/,  
     + RCwf(32,:)/2.803421584636580d0,1.17712d0,1.829324563361060d0/,  
     + RCwf(33,:)/2.420952200569910d0,-0.00000d0,2.594263331494400d0/, 
     + RCwf(34,:)/2.803421584636580d0,-1.17712d0,1.829324563361060d0/, 
     + RCwf(35,:)/1.985812298697740d0,-2.30246d0,1.829324563361060d0/, 
     + RCwf(36,:)/1.749433219688310d0,-3.02996d0,0.591627637208844d0/, 
     + RCwf(37,:)/0.365645988478857d0,-3.47958d0,0.591627637208844d0/, 
     + RCwf(38,:)/-0.365645988478858d0,-3.47958d0,-0.591627637208844d0/
     + RCwf(39,:)/-1.74943321968831d0,-3.02996d0,-0.591627637208844d0/,
     + RCwf(40,:)/-2.34106085689715d0,-2.60012d0,0.591627637208844d0/, 
     + RCwf(41,:)/-3.19628839898270d0,-1.42300d0,0.591627637208844d0/, 
     + RCwf(42,:)/-3.42227004771268d0,-0.72750d0,-0.591627637208844d0/,
     + RCwf(43,:)/-3.42227004771268d0,0.72750d0,-0.591627637208844d0/, 
     + RCwf(44,:)/-3.19628839898270d0,1.42300d0,0.591627637208844d0/,  
     + RCwf(45,:)/-2.34106085689715d0,2.60012d0,0.591627637208844d0/,  
     + RCwf(46,:)/-1.57612208876381d0,2.60012d0,1.829324563361060d0/,  
     + RCwf(47,:)/-0.253202474597252d0,3.02996d0,1.829324563361060d0/, 
     + RCwf(48,:)/0.748115372545527d0,2.30246d0,2.594263331494400d0/,  
     + RCwf(49,:)/0.382469384066670d0,1.17712d0,3.325555308452120d0/,  
     + RCwf(50,:)/1.237696926152220d0,0.00000d0,3.325555308452120d0/,  
     + RCwf(51,:)/0.382469384066670d0,-1.17712d0,3.325555308452120d0/, 
     + RCwf(52,:)/0.748115372545527d0,-2.30246d0,2.594263331494400d0/, 
     + RCwf(53,:)/-0.253202474597252d0,-3.02996d0,1.829324563361060d0/,
     + RCwf(54,:)/-1.576122088763810d0,-2.60012d0,1.829324563361060d0/,
     + RCwf(55,:)/-1.958591472830480d0,-1.42300d0,2.594263331494400d0/,
     + RCwf(56,:)/-2.959909319973260d0,-0.69550d0,1.829324563361060d0/,
     + RCwf(57,:)/-2.959909319973260d0,0.69550d0,1.829324563361060d0/, 
     + RCwf(58,:)/-1.958591472830480d0,1.42300d0,2.594263331494400d0/, 
     + RCwf(59,:)/-1.001317847142780d0,0.72750d0,3.325555308452120d0/, 
     + RCwf(60,:)/-1.001317847142780d0,-0.72750d0,3.325555308452120d0/

c     
c     prepare rotational matrix for water 
c     obtain the SFF coordinates for H1, H2, and O of water 
      call matpre(Eulang, rotmat)
      do i=1,3
         RO_sf(i)=0.d0
      enddo
      call rottrn(rotmat, ROwf, RO_sf, com)
      do i=1,3
         RH1_sf(i)=0.d0
      enddo
      call rottrn(rotmat, RH1wf, RH1_sf, com)
      do i=1,3
         RH2_sf(i)=0.d0
      enddo
      call rottrn(rotmat, RH2wf, RH2_sf, com)
      

      do j=1,60
       do i=1,3
      RC_sf(j,i)=RCwf(j,i)
       enddo
      enddo
       

c writing the coordinates of H2OC60 in SFF
c      do i=1,60
c      write(44,*) RC_sf(i,1:3) 
c      enddo
c      write(44,*) RO_sf(1:3),RH1_sf(1:3), RH2_sf(1:3)     
c 
c ... calculate water dimer energies through LJ potential
      E_H2OC60=0.d0


c ... C-O interaction
      v_colj=0.d0
      do j=1,60  !loop over number of carbon atoms 
      rco=0.d0
        do i=1,3 !loop over x,y,z
        rco=rco+(RO_sf(i)-RC_sf(j,i))*(RO_sf(i)-RC_sf(j,i))
        enddo ! loop over x,y,z
      rco=dsqrt(rco)


      rco6=rco**6.d0
      rco12=rco**12.d0
      v_colj=v_colj+4.d0*epsco*(sigco**12.d0/rco12-sigco**6d0/rco6)
      v=4.d0*epsco*(sigco**12.d0/rco12-sigco**6d0/rco6)
c      write(44,*)'j rco ',j, rco,v
      enddo ! loop over C60


c ... C-H1 interaction
      v_ch1lj=0.d0
      do j=1,60  !loop over number of carbon atoms 
      rch1=0.d0
        do i=1,3 !loop over x,y,z
        rch1=rch1+(RH1_sf(i)-RC_sf(j,i))*(RH1_sf(i)-RC_sf(j,i))
        enddo ! loop over x,y,z
      rch1=dsqrt(rch1)
      rch16=rch1**6.d0
      rch112=rch1**12.d0
      v_ch1lj=v_ch1lj+4.d0*epsch*(sigch**12.d0/rch112-sigch**6d0/rch16)
      v=4.d0*epsch*(sigch**12.d0/rch212-sigch**6d0/rch16)
c      write(44,*)'j rch1 ',j, rch1,v
      enddo ! loop over C60


c ... C-H2 interaction
      v_ch2lj=0.d0
      do j=1,60  !loop over number of carbon atoms 
      rch2=0.d0
        do i=1,3 !loop over x,y,z
        rch2=rch2+(RH2_sf(i)-RC_sf(j,i))*(RH2_sf(i)-RC_sf(j,i))
        enddo ! loop over x,y,z
      rch2=dsqrt(rch2)
      rch26=rch2**6.d0
      rch212=rch2**12.d0
      v=4.d0*epsch*(sigch**12.d0/rch212-sigch**6d0/rch26)
      v_ch2lj=v_ch2lj+4.d0*epsch*(sigch**12.d0/rch212-sigch**6d0/rch26)
c      write(44,*)'j rch2  ',j, rch2, v

      enddo ! loop over C60

      E_H2OC60=(v_colj+v_ch1lj+v_ch2lj)!*kcal2k
      write(44,71) RO_sf(1:3),RH1_sf(1:3), RH2_sf(1:3), E_H2OC60,
     +E_H2OC60/kcal2k  
   71 format(9F12.6,2F18.4)
c      write(44,*)'E_H2OC60=',E_H2OC60

      return
      end


c-----------------------------------------------------------------------
      subroutine matpre(Eulang,rotmat)
      implicit double precision(a-h,o-z)

      dimension Eulang(3),rotmat(3,3)

      phi=Eulang(1)
      theta=Eulang(2)
      chi=Eulang(3)

      cp=cos(phi)
      sp=sin(phi)
      ct=cos(theta)
      st=sin(theta)
      ck=cos(chi)
      sk=sin(chi)

      rotmat(1,1)=cp*ct*ck-sp*sk
      rotmat(1,2)=-cp*ct*sk-sp*ck
      rotmat(1,3)=cp*st
      rotmat(2,1)=sp*ct*ck+cp*sk
      rotmat(2,2)=-sp*ct*sk+cp*ck
      rotmat(2,3)=sp*st
      rotmat(3,1)=-st*ck
      rotmat(3,2)=st*sk
      rotmat(3,3)=ct

      return
      end
c-----------------------------------------------------------------------
      subroutine rottrn(rotmat,rwf,rsf,rcom)
      implicit double precision(a-h,o-z)
      dimension rotmat(3,3),rwf(3),rsf(3),rcom(3)

      do i=1,3
        rsf(i)=rcom(i)
        do j=1,3
          rsf(i)=rsf(i)+rotmat(i,j)*rwf(j)
        enddo
      enddo

      return
      end
c-----------------------------------------------------------------------


