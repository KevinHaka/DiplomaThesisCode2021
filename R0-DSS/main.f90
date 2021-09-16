Program R0_DSS
    implicit none

    !initial conditions and parameters
    integer,   parameter :: NMC=1e2           !Number Monte Carlo
    integer,   parameter :: x=30,y=30         !lattice(x*y=population)
    real,      parameter :: dc=6.             !daily checks (dc->infinity => dt->0)
    character, parameter :: Lp='S'            !shape of interaction, R=Rhombus, C=Circle, S=Square (Lp spaces)
    character, parameter :: PBCs='N'          !Periodic Boundary Conditions, (Y=Yes,N=No)
    character, parameter :: visualization='N' !Ovito visualization results (Y=Yes,N=No)

    character(3), parameter :: variable="inf"     !variable on the x axis (dis=distance, dcp, inf=infectiveness)
    integer,      parameter :: np=21              !number of points (R0) (np>=2))
    real,         parameter :: x_min=0.,x_max=1.  !range of x
    integer,      parameter :: no=2               !which point of the curve to be visualized (1<=no<=np)

    integer :: distance=2         !movement distance
    real    :: dcp=5.             !daily contacts per person
    real    :: infectiveness=0.2  !possibility of infection per contact

    !Áttention, the probability of stay is calculated by the function P=t/(1+t)
    real, parameter :: C(3)=(/1.,0.,0./)      !stochastic vector (C=Cured)
    real, parameter :: I(3)=(/0.2,0.75,0.05/) !stochastic vector (I=Infected)
    real, parameter :: D(3)=(/0.,0.,1./)      !stochastic vector (R=Removed)

    !do not change these definitions
    integer, parameter :: inps(3)=(/0,1,0/)        !initial number of people in each situation
    character          :: names(3)=(/'C','I','D'/) !names of situations
    character          :: space(x,y)               !space of individuals and their condition
    character          :: cp(2)                    !condition of the person 1,2
    integer            :: pp(2,2)                  !position of person 1,2 (1=x,2=y)
    integer            :: alive                   !alive persons
    integer            :: contacts                 !contacts per check
    integer            :: ci,cj,ck,cq,cw,cr        !counters
    integer            :: random_integer           !function
    integer            :: targ(2)                  !target (1=x,2=y)
    real               :: ra_num                   !random number
    real               :: corrector                !rounding corrector
    real               :: Sm(3,3)                  !Stochastic matrix
    real               :: R0(NMC)                  !basic reproduction number
    real               :: av_R0                    !average basic reproduction number
    real               :: e_R0                     !error of R0

    call init_random_seed

    open(1,file="Ovito.dat")
    open(2,file="Results.dat")

    !creating a stochastic matrix
    Sm(1,:)=C
    Sm(2,:)=I
    Sm(3,:)=D

    !transformation of probabilities
    do cq=1,size(names)
        if (SM(cq,cq)==1.) go to 1
        SM(cq,cq)=dc*SM(cq,cq)/(1-SM(cq,cq)+dc*SM(cq,cq))
        ra_num=SM(cq,cq)
        sm(cq,:)=(1-SM(cq,cq))*sm(cq,:)/(sum(SM(cq,:))-SM(cq,cq))
        SM(cq,cq)=ra_num
1   end do

    !calculations of R0
    do cr=0,np-1
        R0=0.

        !x-axis
        selectcase (variable)
        case ("dis")
            distance=cr*nint(x_max-x_min)/(np-1)+nint(x_min)
            print*, distance
        case ("dcp")
            dcp=cr*(x_max-x_min)/(np-1)+x_min
            print*, dcp
        case ("inf")
            infectiveness=cr*(x_max-x_min)/(np-1)+x_min
            print*, infectiveness
        end select

        !Monte Carlo calculations
        do cq=1,NMC
            space='H'
            alive=x*y
            corrector=0.

            !random infections in the population
            do cw=1,size(names)
                do ck=1,inps(cw)
2                   pp(1,:)=(/random_integer(1,x),random_integer(1,y)/)

                    if (space(pp(1,1),pp(1,2))/='H') go to 2

                    space(pp(1,1),pp(1,2))=names(cw)

                    if (cw==2) targ(:)=pp(1,:)
                end do
            end do

            !Ovito visualization
            if ((cq==1).and.(cr==no).and.(visualization=='Y')) call Ovito(x,y,space,len(space),1)

            !running time
            do
                !contacts number calculation
                contacts=nint(dcp*alive/(2*dc)+corrector)
                corrector=corrector+dcp*alive/(2*dc)-contacts

                !random contacts and infections
                do cw=1,contacts
                    call random_contact(x,y,distance,Lp,pp,cp,space,PBCs)

                    !ban meeting with dead
                    if ((cp(1)=='D').or.(cp(2)=='D')) corrector=corrector+1

                    !healthy-infected contact
                    if ((cp(1)=='H').and.(cp(2)=='I')) then
                        call random_number(ra_num)

                        if (ra_num<infectiveness) then
                            space(pp(1,1),pp(1,2))='I'

                            if (all(targ(:)==pp(2,:))) R0(cq)=R0(cq)+1.
                        endif

                        elseif ((cp(1)=='I').and.(cp(2)=='H')) then
                            call random_number(ra_num)

                            if (ra_num<infectiveness) then
                                space(pp(2,1),pp(2,2))='I'

                                if (all(targ(:)==pp(1,:))) R0(cq)=R0(cq)+1.
                            endif
                    end if
                end do

                !disease progression in any infected
                do ci=1,x
                    do cj=1,y
                        call random_number(ra_num)

                        selectcase (space(ci,cj))
                            case ('I')
                                do cw=1,3
                                    if (ra_num<Sum(Sm(2,1:cw))) then
                                        space(ci,cj)=names(cw)

                                        if (names(cw)=='D') alive=alive-1
                                        exit
                                    endif
                                end do
                        end select
                    end do
                end do

                !Ovito visualization
                if ((cq==1).and.(cr==no).and.(visualization=='Y')) call Ovito(x,y,space,len(space),1)

                if ((space(targ(1),targ(2))/="I")) exit
                if (alive<=1) exit
            end do
        end do

        av_R0=sum(R0)/NMC
        e_R0=sqrt(sum((R0-av_R0)**2)/NMC/(NMC-1))

        !significant digits
        if (e_R0/=0.) then
            ck=1

            do while (e_R0<=1.)
                e_R0=e_R0*10
                ck=ck*10
            end do

            e_R0=real(nint(e_R0))/ck
            if (e_R0==10./ck) ck=ck/10
            av_R0=real(nint(av_R0*ck))/ck
        end if

        !x-axis
        selectcase (variable)
        case ("dis")
            write(2,*) distance,av_R0,e_R0
        case ("dcp")
            write(2,*) dcp,av_R0,e_R0
        case ("inf")
            write(2,*) infectiveness,av_R0,e_R0
        end select
    end do
end program R0_DSS


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
    integer,   intent(in)  :: x,y            !dimension lengths
    integer,   intent(in)  :: distance       !movement distance
    integer,   intent(out) :: pp(2,2)        !position of person 1,2 (1=x,2=y)
    character, intent(out) :: cp(2)          !condition of the person 1,2
    character, intent(in)  :: space(x,y)     !space of individuals and their condition
    character, intent(in)  :: Lp             !shape of interaction, R=Rhombus, C=Circle, S=Square (Lp spaces)
    character, intent(in)  :: PBCs           !periodic boundary conditions
    integer                :: random_integer !function
    real                   :: n              !variable

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
