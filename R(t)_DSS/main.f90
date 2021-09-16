Program R_DSS
    implicit none

    !initial conditions and parameters
    integer,   parameter :: NMC=5e1           !Number Monte Carlo (NMC/=1)
    integer,   parameter :: x=50,y=50         !lattice(x*y=population)
    integer,   parameter :: of=12*2*30        !Ovito frames (of<day*dc)
    integer,   parameter :: days=2*30         !time
    real,      parameter :: dc=12.            !daily checks (dc->infinity => dt->0)
    character, parameter :: Lp='S'            !shape of interaction, R=Rhombus, C=Circle, S=Square (Lp spaces)
    character, parameter :: PBCs='Y'          !Periodic Boundary Conditions, (Y=Yes,N=No)
    character, parameter :: visualization='N' !Ovito visualization results (Y=Yes,N=No)

    character(3), parameter :: variable="inf"
    !parameter change (inc=initial number of cured, dis=distance, dcp, inf=infectiveness)
    integer, parameter      :: np=2                !number of curves R(t) (np=1,2,3,4)
    real,    parameter      :: x_min=0.3,x_max=0.6 !parameter values

    integer :: inps(5)=(/0,1,0,0,0/) !initial number of people in each situation
    integer :: distance=2            !movement distance
    real    :: dcp=4.                !daily contacts per person
    real    :: infectiveness=0.3     !possibility of infection per contact

    !Áttention, the probability of stay is calculated by the function P=t/(1+t)
    character(len=2) :: names(5)=(/'In','Im','Is','C ','D '/) !names of situations
    real,  parameter :: In(5)=(/0.5,0.5,0.,0.,0./)            !stochastic vector (In=Incubation)
    real,  parameter :: Im(5)=(/0.,2/3.,1/12.,1/4.,0./)       !stochastic vector (Im=Infected mild)
    real,  parameter :: Is(5)=(/0.,0.,0.75,0.2,0.05/)         !stochastic vector (Is=Infected severe)
    real,  parameter :: C(5)=(/0.,0.,0.,1.,0./)               !stochastic vector (C=Cured)
    real,  parameter :: D(5)=(/0.,0.,0.,0.,1./)               !stochastic vector (D=Dead)

    !do not change these definitions
    integer, parameter :: tc=floor(days*dc)        !total checks
    character(2)       :: space(x,y)               !space of individuals and their condition
    character(2)       :: cp(2)                    !condition of the person 1,2
    integer            :: pp(2,2)                  !position of person 1,2 (1=x,2=y)
    integer            :: contacts                 !contacts per check
    integer            :: ci,cj,ck,cq,cw,cr        !counters
    integer            :: random_integer           !function
    real               :: ra_num                   !random number
    real               :: corrector                !rounding corrector
    real               :: Sm(5,5)                  !Stochastic matrix
    real               :: R(0:tc,NMC)              !effective reproduction number
    real               :: av_R(0:tc)               !average effective reproduction number
    real               :: e_R(0:tc)                !error of R
    integer            :: helper(x,y,2)            !auxiliary table
    !1=number of people he infects, 2=time at which become infected
    integer            :: helper2(0:tc,2)          !auxiliary table 2
    !1=number of people infected, 2=number of people infecting
    integer            :: freq(0:tc,NMC)           !frequency of at least one infection in a given day
    real               :: av_freq(0:tc)            !average frequency of at least one infection in a given day
    real               :: e_freq(0:tc)             !error of frequency of at least one infection in a given day
    integer            :: Npes(0:tc,6,NMC)         !number of people in each situation at any time (1=In, 2=Im, 3=Is, 4=C, 5=D, 6=H)
    real               :: av_Npes(0:tc,6)          !average number of people in each situation at any time
    real               :: e_Npes(0:tc,6)           !error of number of people in each situation at any time
    integer            :: Ped(0:tc,NMC)            !Probability of existing a disease at a particular day
    real               :: av_Ped(0:tc)             !average Probability of existing a disease at a particular day
    real               :: e_Ped(0:tc)              !error of Probability of existing a disease at a particular day
    logical            :: TF                       !controller

    call init_random_seed

    !Ovito data
    open(41,file="Ovito1.dat")
    open(42,file="Ovito2.dat")
    open(43,file="Ovito3.dat")
    open(44,file="Ovito4.dat")

    !av_R, e_R
    open(1,file="R1.dat")
    open(2,file="R2.dat")
    open(3,file="R3.dat")
    open(4,file="R4.dat")

    !av_Npes, e_Npes
    open(11,file="R11.dat")
    open(12,file="R12.dat")
    open(13,file="R13.dat")
    open(14,file="R14.dat")

    !av_Ped, e_Ped
    open(21,file="R21.dat")
    open(22,file="R22.dat")
    open(23,file="R23.dat")
    open(24,file="R24.dat")

    !av_freq, e_freq
    open(31,file="R31.dat")
    open(32,file="R32.dat")
    open(33,file="R33.dat")
    open(34,file="R34.dat")

    !creating a stochastic matrix
    Sm(1,:)=In
    Sm(2,:)=Im
    Sm(3,:)=Is
    Sm(4,:)=C
    Sm(5,:)=D

    !transformation of probabilities
    do cq=1,size(names)
        if (Sm(cq,cq)==1.) go to 1
        Sm(cq,cq)=dc*Sm(cq,cq)/(1+Sm(cq,cq)*(dc-1.))
        ra_num=Sm(cq,cq)
        sm(cq,:)=(1-Sm(cq,cq))*Sm(cq,:)/(sum(Sm(cq,:))-Sm(cq,cq))
        Sm(cq,cq)=ra_num
1   end do

    !calculations of R(t)
    do cr=1,np
        R=0.
        Npes=0
        Ped=0
        freq=0

        !parameter change
        selectcase (variable)
        case ("dis")
            if (np>=2) then
                distance=(cr-1)*nint(x_max-x_min)/(np-1)+nint(x_min)

                else
                    distance=nint(x_min)
            end if
            print*, cr,distance

        case ("dcp")
            if (np>=2) then
                dcp=(cr-1)*(x_max-x_min)/(np-1)+x_min

                else
                    dcp=x_min
            end if
            print*, cr,dcp

        case ("inf")
            if (np>=2) then
                infectiveness=(cr-1)*(x_max-x_min)/(np-1)+x_min

                else
                    infectiveness=x_min
            end if
            print*, cr,infectiveness

        case ("inc")
            if (np>=2) then
                inps(4)=(cr-1)*nint(x_max-x_min)/(np-1)+nint(x_min)

                else
                    inps(4)=nint(x_min)
            end if
            print*, cr,inps
        end select

        !Monte Carlo calculations
        do cq=1,NMC
            space='H'
            corrector=0.
            helper(:,:,1)=0
            helper(:,:,2)=-1
            helper2=0
            Npes(0,:,cq)=(/inps(:),x*y-sum(inps(:))/)
            Ped(0,cq)=1
            freq(0,cq)=1

            !random infections in the population
            do cw=1,size(names)
                do ck=1,inps(cw)
2                   pp(1,:)=(/random_integer(1,x),random_integer(1,y)/)

                    if (space(pp(1,1),pp(1,2))/='H') go to 2

                    space(pp(1,1),pp(1,2))=names(cw)

                    if (cw<=3) helper(pp(1,1),pp(1,2),2)=0
                end do
            end do

            !Ovito visualization
            if ((cq==1).and.(visualization=='Y')) call Ovito(x,y,space,len(space),cr+40)

            !running time
            do ck=1,tc
                TF=.True.
                Npes(ck,:,cq)=Npes(ck-1,:,cq)
                if (sum(Npes(ck,1:3,cq))>=1) Ped(ck,cq)=Ped(ck,cq)+1

                !contacts number calculation
                contacts=nint(dcp*(x*y-inps(5))/(2*dc)+corrector)
                corrector=corrector+dcp*x*(y-inps(5))/(2*dc)-contacts

                !random contacts and infections
                do cw=1,contacts
                    if (x*y-Npes(ck,5,cq)<=1) exit
3                   call random_contact(x,y,distance,Lp,pp,cp,space,PBCs)

                    !ban meeting with dead
                    if ((cp(1)=='D').or.(cp(2)=='D')) go to 3

                    !healthy-infected contact
                    if ((cp(1)=='H').and.((cp(2)=='Im').or.(cp(2)=='Is'))) then
                        call random_number(ra_num)

                        if (ra_num<infectiveness) then
                            space(pp(1,1),pp(1,2))='In'
                            Npes(ck,1,cq)=Npes(ck,1,cq)+1
                            Npes(ck,6,cq)=Npes(ck,6,cq)-1

                            if (TF) then
                                freq(ck,cq)=freq(ck,cq)+1
                                TF=.False.
                            end if

                            helper(pp(1,1),pp(1,2),2)=ck
                            helper(pp(2,1),pp(2,2),1)=helper(pp(2,1),pp(2,2),1)+1
                        endif

                        elseif (((cp(1)=='Im').or.(cp(1)=='Is')).and.(cp(2)=='H')) then
                            call random_number(ra_num)

                            if (ra_num<infectiveness) then
                                space(pp(2,1),pp(2,2))='In'
                                Npes(ck,1,cq)=Npes(ck,1,cq)+1
                                Npes(ck,6,cq)=Npes(ck,6,cq)-1

                                if (TF) then
                                    freq(ck,cq)=freq(ck,cq)+1
                                    TF=.False.
                                end if

                                helper(pp(2,1),pp(2,2),2)=ck
                                helper(pp(1,1),pp(1,2),1)=helper(pp(1,1),pp(1,2),1)+1
                            endif
                    end if
                end do

                !disease progression in any infected
                do ci=1,x
                    do cj=1,y
                        call random_number(ra_num)

                        selectcase (space(ci,cj))
                            case ('In')
                                do cw=1,size(names)
                                    if (ra_num<Sum(SM(1,1:cw))) then
                                        space(ci,cj)=names(cw)
                                        Npes(ck,cw,cq)=Npes(ck,cw,cq)+1
                                        Npes(ck,1,cq)=Npes(ck,1,cq)-1
                                        exit
                                    endif
                                end do

                            case ('Im')
                                do cw=1,size(names)
                                    if (ra_num<Sum(SM(2,1:cw))) then
                                        space(ci,cj)=names(cw)
                                        Npes(ck,cw,cq)=Npes(ck,cw,cq)+1
                                        Npes(ck,2,cq)=Npes(ck,2,cq)-1
                                        exit
                                    endif
                                end do

                            case ('Is')
                                do cw=1,size(names)
                                    if (ra_num<Sum(SM(3,1:cw))) then
                                        space(ci,cj)=names(cw)
                                        Npes(ck,cw,cq)=Npes(ck,cw,cq)+1
                                        Npes(ck,3,cq)=Npes(ck,3,cq)-1
                                        exit
                                    endif
                                end do

                            case ('C')
                                do cw=1,size(names)
                                    if (ra_num<Sum(SM(4,1:cw))) then
                                        space(ci,cj)=names(cw)
                                        Npes(ck,cw,cq)=Npes(ck,cw,cq)+1
                                        Npes(ck,4,cq)=Npes(ck,4,cq)-1
                                        exit
                                    endif
                                end do

                            case ('D')
                                do cw=1,size(names)
                                    if (ra_num<Sum(SM(5,1:cw))) then
                                        space(ci,cj)=names(cw)
                                        Npes(ck,cw,cq)=Npes(ck,cw,cq)+1
                                        Npes(ck,5,cq)=Npes(ck,5,cq)-1
                                        exit
                                    endif
                                end do
                        end select
                    end do
                end do

                !Ovito visualization
                if ((cq==1).and.(visualization=='Y').and.(mod(ck,tc/of)==0)) call Ovito(x,y,space,len(space),cr+40)
            end do

            !calculations
            do ci=1,x
                do cj=1,x
                    if (helper(ci,cj,2)>=0) then
                        helper2(helper(ci,cj,2),1)=helper2(helper(ci,cj,2),1)+helper(ci,cj,1)
                        helper2(helper(ci,cj,2),2)=helper2(helper(ci,cj,2),2)+1
                    end if
                end do
            end do

            do ck=0,tc
                if (helper2(ck,2)/=0) R(ck,cq)=helper2(ck,1)/real(helper2(ck,2))
            end do
        end do

        !calculation of mean values ​​and errors
        do cw=0,tc
            av_Ped(cw)=sum(Ped(cw,:))/real(NMC)
            e_Ped(cw)=sqrt(sum((Ped(cw,:)-av_Ped(cw))**2)/NMC/(NMC-1))

            av_freq(cw)=sum(freq(cw,:))/real(NMC)
            e_freq(cw)=sqrt(sum((freq(cw,:)-av_freq(cw))**2)/NMC/(NMC-1))

            do ck=1,6
                av_Npes(cw,ck)=sum(Npes(cw,ck,:))/real(NMC)
                e_Npes(cw,ck)=sqrt(sum((Npes(cw,ck,:)-av_Npes(cw,ck))**2)/NMC/(NMC-1))
            end do

            if (av_freq(cw)==0.) then
                av_R(cw)=0.
                e_R(cw)=0.

                else if (av_freq(cw)*NMC<=5.) then
                    av_R(cw)=sum(R(cw,:))/av_freq(cw)/NMC
                    e_R(cw)=0.

                    else
                        av_R(cw)=sum(R(cw,:))/av_freq(cw)/NMC
                        e_R(cw)=sqrt(sum((R(cw,:)-av_R(cw))**2)/(av_freq(cw)*NMC)/((av_freq(cw)*NMC)-1.))
            end if
        end do

        !results, per day point
        write(cr,*) 0,av_R(0),e_R(0)
        write(cr+20,*) 0,av_Ped(0),e_Ped(0)
        write(cr+30,*) 0,av_freq(0),e_freq(0)
        do ck=0,days-1
            write(cr,*) ck+1,sum(av_R(ck*nint(dc)+1:(ck+1)*nint(dc)))/dc,sqrt(sum(e_R(ck*nint(dc)+1:(ck+1)*nint(dc))**2))/dc
            write(cr+20,*) ck+1,sum(av_Ped(ck*nint(dc)+1:(ck+1)*nint(dc)))/dc,sqrt(sum(e_Ped(ck*nint(dc)+1:(ck+1)*nint(dc))**2))/dc
            write(cr+30,*)ck+1,sum(av_freq(ck*nint(dc)+1:(ck+1)*nint(dc)))/dc,&
            sqrt(sum(e_freq(ck*nint(dc)+1:(ck+1)*nint(dc))**2))
        end do

        !results
        do ck=0,tc
            !write(cr,*) ck,av_R(ck),e_R(ck)
            write(cr+10,*) ck,av_Npes(ck,:),e_Npes(ck,:)
            !write(cr+20,*) ck,av_Ped(ck),e_Ped(ck)
            !write(cr+30,*) ck,av_freq(ck),e_freq(ck)
        end do
    end do
end program R_DSS


integer function random_integer(a,b)
    !random integer value
    implicit none
    real                :: ra_var !random variable
    integer, intent(in) :: a,b    !space limits

    call random_number(ra_var)
    random_integer=floor((b-a+1)*ra_var)+a
end function


subroutine random_contact(x,y,distance,Lp,pp,cp,space,PBCs)
    !random allowed movements and contacts
    implicit none
    integer,      intent(in)  :: x,y            !dimension lengths
    integer,      intent(in)  :: distance       !movement distance
    integer,      intent(out) :: pp(2,2)        !position of person 1,2 (1=x,2=y)
    character(2), intent(out) :: cp(2)          !condition of the person 1,2
    character(2), intent(in)  :: space(x,y)     !space of individuals and their condition
    character,    intent(in)  :: Lp             !shape of interaction, R=Rhombus, C=Circle, S=Square (Lp spaces)
    character,    intent(in)  :: PBCs           !periodic boundary conditions
    integer                   :: random_integer !function
    real                      :: n              !variable

    selectcase (Lp)
    case ('R')
        n=1.
    case ('C')
        n=2.
    endselect

    !first random person
    pp(1,1)=random_integer(1,x)
    pp(1,2)=random_integer(1,y)
    cp(1)=space(pp(1,1),pp(1,2))

    !periodic boundary conditions
    selectcase (PBCs)
    case ('Y')
        !second random person at allowed distance
        do
            pp(2,1)=random_integer(pp(1,1)-distance,pp(1,1)+distance)
            pp(2,2)=random_integer(pp(1,2)-distance,pp(1,2)+distance)

            if ((Lp=='S').and.(abs(pp(2,1)-pp(1,1))+abs(pp(2,2)-pp(1,2))/=0)) go to 2

            if ((distance**n>=abs(pp(2,1)-pp(1,1))**n+abs(pp(2,2)-pp(1,2))**n) &
            .and.(abs(pp(2,1)-pp(1,1))+abs(pp(2,2)-pp(1,2))/=0)) then

                !for the x
2               if (pp(2,1)<1) then
                    pp(2,1)=mod(pp(2,1),x)
                    pp(2,1)=pp(2,1)+x

                    elseif (pp(2,1)>x) then
                        pp(2,1)=mod(pp(2,1),x)

                        if (pp(2,1)==0) pp(2,1)=x
                end if

                !for the y
                if (pp(2,2)<1) then
                    pp(2,2)=mod(pp(2,2),y)
                    pp(2,2)=pp(2,2)+y

                    elseif (pp(2,2)>y) then
                        pp(2,2)=mod(pp(2,2),y)

                        if (pp(2,2)==0) pp(2,2)=y
                end if

                cp(2)=space(pp(2,1),pp(2,2))
                exit
            end if
        end do

    case ('N')
        !second random person at allowed distance
        do
            !for the x
            if ((pp(1,1)-distance<1).and.(pp(1,1)+distance>x)) then
                pp(2,1)=random_integer(1,x)

                elseif (pp(1,1)-distance<1) then
                    pp(2,1)=random_integer(1,pp(1,1)+distance)

                elseif (pp(1,1)+distance>x) then
                    pp(2,1)=random_integer(pp(1,1)-distance,x)

                    else
                        pp(2,1)=random_integer(pp(1,1)-distance,pp(1,1)+distance)
            end if

            !for the y
            if ((pp(1,2)-distance<1).and.(pp(1,2)+distance>y)) then
                pp(2,2)=random_integer(1,y)

                elseif (pp(1,2)-distance<1) then
                    pp(2,2)=random_integer(1,pp(1,2)+distance)

                elseif (pp(1,2)+distance>y) then
                    pp(2,2)=random_integer(pp(1,2)-distance,y)

                    else
                        pp(2,2)=random_integer(pp(1,2)-distance,pp(1,2)+distance)
            end if

            if ((Lp=='S').and.(abs(pp(2,1)-pp(1,1))+abs(pp(2,2)-pp(1,2))/=0)) goto 3

            if ((distance**n>=abs(pp(2,1)-pp(1,1))**n+abs(pp(2,2)-pp(1,2))**n) &
            .and.(abs(pp(2,1)-pp(1,1))+abs(pp(2,2)-pp(1,2))/=0)) then

3               cp(2)=space(pp(2,1),pp(2,2))
                exit
            end if
        end do
    endselect
end subroutine


subroutine Ovito(x,y,space,nsl,nn)
    !Ovito visualization
    implicit none
    integer,            intent(in) :: nn         !file number (open in the main program)
    integer,            intent(in) :: nsl        !number of space len
    integer,            intent(in) :: x,y        !dimension lengths
    character(len=nsl), intent(in) :: space(x,y) !space of individuals and their condition
    integer                        :: i,j        !counters

    write(nn,*) x*y
    write(nn,*)
    do i=1,x
        do j=1,y
            write(nn,*) space(i,j),i,j
        end do
    end do
end subroutine


subroutine init_random_seed()
! random_seed based on system time
    implicit none
    integer :: i, n, clock
    integer, dimension(:), allocatable :: seed

    call random_seed(size = n)
    allocate(seed(n))

    call system_clock(count=clock)

    seed = clock + 37 * (/ (i - 1, i = 1, n) /)
    call random_seed(put = seed)

    deallocate(seed)
end subroutine init_random_seed

