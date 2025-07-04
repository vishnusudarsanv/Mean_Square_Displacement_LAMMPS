program msd

!trajectory reader parameter variables
implicit none
real(kind=8),allocatable::na(:,:,:),fsi(:,:,:)
real::xlo,xhi,ylo,yhi,zlo,zhi,boxx,boxy,boxz,x,y,z,volume,boxzz
real:: sx,sy,sz
integer::junk_num,ok,totbinnumd,omm,eq_pos,ios
integer::nframes,natoms,i,j,atmid,atmty,atomno,mol,mol_na,mol_fsi,fcount
character(len=200) :: line
character*200::trj_filename,system_path,value_str,trj_filepath,output_filepath,output_file
character*10::junk_char
real, parameter :: PI=4.D0*DATAN(1.D0)
integer, parameter ::totframes=4196,nlip=110      ! 450 for LAGP

!histogram parameter variables
real:: normFac, rr2, dx, dy, dz, rr,junk2,delR2,rangenumd,dist_int,r_iz,l_iz,temp_na
real,allocatable::msd_x(:),msd_x_part(:),msd_y(:),msd_y_part(:),msd_z(:),msd_z_part(:),msd_tot(:),msd_xy(:),msd_xy_part(:),msd_tot_accum(:)
real,allocatable:: msd_corr(:),msd_corr_part(:,:),norm(:),surv(:)
integer:: sod, k, n, atyp1, totbin, nframe,junk1,sol,nfsi,l,q,h,nbin,Nik,k1,k2,m,temp_ref
integer, allocatable:: ref_id(:,:),na_nbin(:),counter(:)

integer:: maxtimeorigin,tmax,nsod,lip,lig,Ntt_sum
integer, allocatable:: li_id(:)
!this section is a trajectory file reader

!------------------------Config parameters -----------------!
!--------Config file parameters-------!

    system_path = ""
    trj_filename =""
    dist_int = 0
    l_iz = 0
    r_iz=0.0
!-------------------------------------!
open(unit=11, file="inp_msd_tot.dat")
!-----------------------------------------!
    do 20
        read(11, '(A)', iostat=ok) line
        if (ok<0) exit

        if (index(line, 'system_path') > 0) then
            eq_pos = index(line, '=')
            if (eq_pos > 0) then
                system_path = adjustl(trim(line(eq_pos+1:)))
            endif
        endif

        if (index(line, 'trj_filename') > 0) then
            eq_pos = index(line, '=')
            if (eq_pos > 0) then
                trj_filename = adjustl(trim(line(eq_pos+1:)))
                trj_filename = trj_filename(2:len_trim(trj_filename)-1)
            endif
        endif
        if (index(line, 'output_file') > 0) then
            eq_pos = index(line, '=')
            if (eq_pos > 0) then
                output_file = adjustl(trim(line(eq_pos+1:)))
                output_file = output_file(2:len_trim(output_file)-1)
            endif
        endif

        if (index(line, 'lip') > 0) then
            eq_pos = index(line, '=')
            if (eq_pos > 0) then
                value_str = trim(line(eq_pos+1:))
                read(value_str, *, iostat=ok) lip
                if (ok /= 0) lip = 0
            endif
        endif
        if (index(line, 'lig') > 0) then
            eq_pos = index(line, '=')
            if (eq_pos > 0) then
                value_str = trim(line(eq_pos+1:))
                read(value_str, *, iostat=ok) lig
                if (ok /= 0) lig = 0
            endif
        endif        

        if (index(line, 'dist_int') > 0) then
            eq_pos = index(line, '=')
            if (eq_pos > 0) then
                value_str = trim(line(eq_pos+1:))
                read(value_str, *, iostat=ok) dist_int
                if (ok /= 0) dist_int = 0
            endif
        endif

        if (index(line, 'r_iz') > 0) then
            eq_pos = index(line, '=')
            if (eq_pos > 0) then
                value_str = trim(line(eq_pos+1:))
                read(value_str, *, iostat=ok) r_iz
                if (ok /= 0) r_iz = 0
            endif
        endif

        if (index(line, 'l_iz') > 0) then
            eq_pos = index(line, '=')
            if (eq_pos > 0) then
                value_str = trim(line(eq_pos+1:))
                read(value_str, *, iostat=ok) l_iz
                if (ok /= 0) l_iz = 0
            endif
        endif
  20  end do
  close(11)
!----------------------------------------!
trj_filepath=trim(system_path)//trim(trj_filename)
print*, "System Path: ", trim(system_path)
print*, "Trajectory Filename: ", trim(trj_filename)
print*,"Full path",trim(trj_filepath)
print*, "Lipid: ", lip
print*, "Lig : ", lig
print*, "Distance Interval: ", dist_int
print*, "r_iz: ", r_iz
print*, "l_iz: ", l_iz
!--------------------------------

print*,totframes
print*,trj_filepath
output_filepath=trim(system_path)//trim(output_file)
print*,output_filepath

!-----------------------------------------------------------!

!trj_filename=trim("production_comb.lammpstrj")
print*,trj_filename
open(unit=12,file=trj_filepath)
open(unit=15,file=output_filepath)
!open(unit=16,file='diff_20_away_int.dat')
!open(unit=16,file='msd_fsi.dat')
!open(unit=17,file='msd_corr.dat')
read(12,*)junk_char
read(12,*)junk_num
read(12,*)junk_char
read(12,*)natoms
read(12,*)junk_char
read(12,*)xlo,xhi
read(12,*)ylo,yhi
read(12,*)zlo,zhi
read(12,*)junk_char

allocate(na(totframes,3,nlip))
!allocate(ref_id(nlip,totframes))
print*,'number of atoms is',natoms
!close(12)

boxzz=zhi-zlo
tmax= 3700
allocate(ref_id(nlip,totframes))
allocate(msd_x(tmax),msd_x_part(tmax))
allocate(msd_y(tmax),msd_y_part(tmax))
allocate(msd_xy(tmax), msd_xy_part(tmax))
allocate(msd_z(tmax),msd_z_part(tmax))
allocate(msd_tot(tmax), msd_tot_accum(tmax))
allocate(na_nbin(totframes))
allocate(counter(tmax))
allocate(norm(tmax))
allocate(surv(tmax))
allocate(li_id(natoms))
!nframes=0
!open(unit=12,file=trj_filepath)

na=0.0d0
ref_id=0.0d0
rangenumd=boxzz
totbinnumd=500
delR2= rangenumd/(totbinnumd)
na_nbin=0
li_id = 0
k1=0
print*,k1
do i=1,natoms
        read(12,*)atmid,atmty,x,y,z
        if(atmty==lip)then                            !!!!!! selected diffusive atom 
                k1=k1+1
                li_id(atmid)=k1
        endif

enddo
close(12)

nframes=0


open(unit=13,file=trj_filepath)


print*,"Reading frames"
do 31 fcount=1,totframes
    read(13,'(A10)',iostat=ok)junk_char
    if(ok<0)then
            print*,"frames read=",nframes
            exit
    endif
    nframes=nframes+1
    k=0
    read(13,'(I10)')junk_num
    read(13,*)junk_char
    read(13,*)natoms
    read(13,*)junk_char
    read(13,*)xlo,xhi
    read(13,*)ylo,yhi
    read(13,*)zlo,zhi
    read(13,*)junk_char 
    do 32 atomno=1,natoms
        read(13,*)atmid,atmty,x,y,z
        if(atmty==lip)then                               !!!!!! selected diffusive atom 
            na(nframes,1,li_id(atmid)) = x 
            na(nframes,2,li_id(atmid)) = y 
            na(nframes,3,li_id(atmid)) = z 
        end if
        
   32 end do
31 end do

if(nframes<tmax)then
        print*,"Nframes is less than Tmax of diffusion, stopping",nframes,tmax
        stop
endif
!do 21 i=1,nframes
!print*,ref_id(1,i)
!21 end do

k=0
i=0
do 35 i=1,tmax
!print*,'i =',i
        counter(i)=0
        msd_x_part(i)=0.0d0
        msd_y_part(i)=0.0d0
        msd_z_part(i)=0.0d0
        msd_xy_part(i)=0.0d0
        msd_x(i)=0.0d0
        msd_y(i)=0.0d0
        msd_z(i)=0.0d0
        msd_xy(i)=0.0d0
        msd_tot(i)=0.0d0
        msd_tot_accum(i) =0.0d0
        surv(i)=0.0d0
        maxtimeorigin=totframes -i 
        do 36 k1=1,nlip
        !print*,'k1= ',k1
            do 37 k=1,maxtimeorigin,1
            !print*,'k =', k
                    msd_x_part(i)=msd_x_part(i)+(na(k+i,1,k1)-na(k,1,k1))**2
                    !print*,na(k+i,1,k1),na(k,1,k1)
                    !print*, i, 'msd_x_part(i) =',msd_x_part(i)
                    msd_y_part(i)=msd_y_part(i)+(na(k+i,2,k1)-na(k,2,k1))**2
                    msd_z_part(i)=msd_z_part(i)+(na(k+i,3,k1)-na(k,3,k1))**2
                    msd_xy_part(i)=msd_xy_part(i)+0.5*((na(k+i,1,k1)-na(k,1,k1))**2+(na(k+i,2,k1)-na(k,2,k1))**2)
                    msd_tot_accum(i)=msd_tot_accum(i)+(na(k+i,1,k1)-na(k,1,k1))**2+(na(k+i,2,k1)-na(k,2,k1))**2+(na(k+i,3,k1)-na(k,3,k1))**2
                    counter(i)=counter(i)+1
                    !print*,'counter = ',counter(i)

               37  end do
        36 end do
        msd_x(i)=msd_x(i)+msd_x_part(i)
        !print*,'msd_x(i) = ',msd_x(i)
        msd_y(i)=msd_y(i)+msd_y_part(i)
        msd_z(i)=msd_z(i)+msd_z_part(i)
        msd_xy(i)=msd_xy(i)+msd_xy_part(i)
        msd_tot(i)=msd_tot(i)+msd_tot_accum(i)
        norm(i)= real(counter(i))

        write(15,'(2X,F10.7,4X,F15.5,4X,F15.5,4X,F15.5,4X,F15.5,6X,F15.5,6X,F15.5)')i*0.01,msd_x(i)/Dble(counter(i)),msd_y(i)/Dble(counter(i)),msd_z(i)/Dble(counter(i)),msd_xy(i)/Dble(counter(i)),msd_tot(i)/Dble(counter(i)),surv(i)/Dble(maxtimeorigin)
35 end do
!        write(15,'(2X,F10.7,4X,F15.5,4X,F15.5,4X,F15.5,4X,F15.5,6X,F15.5,6X,F15.5)')i*0.01,msd_x(i)/Dble(norm(i)),msd_y(i)/Dble(norm(i)),msd_z(i)/Dble(norm(i)),msd_xy(i)/Dble(norm(i)),msd_tot(i)/Dble(norm(i)),surv(i)/Dble(maxtimeorigin)
!35 end do
close(13)
end program msd
