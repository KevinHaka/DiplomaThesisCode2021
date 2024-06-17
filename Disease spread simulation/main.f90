Program Disease_spread_simulation
    implicit none

    !initial conditions and parameters
    integer,   parameter :: NMC=1e2             !Number Monte Carlo
    integer,   parameter :: x=50,y=50           !lattice(x*y=population)
    integer,   parameter :: distance=1          !movement distance
    integer,   parameter :: of=150              !Ovito frames (of<day*dc)
    integer,   parameter :: days=10             !time
    real,      parameter :: dc=30.              !daily checks (dc->infinity => dt->0)
    real,      parameter :: dcp=4.              !daily contacts per person
    real,      parameter :: infectiveness=0.2   !possibility of infection per contact
    character, parameter :: Lp='S'              !shape of interaction, R=Rhombus, C=Circle, S=Square (Lp spaces)
    character, parameter :: PBCs='Y'            !Periodic Boundary Conditions, (Y=Yes,N=No)
    character, parameter :: visualization='N'   !Ovito visualization results (Y=Yes,N=No)

    !Áttention, the probability of stay is calculated by the function P=t/(1+t)
    character(len=2)   :: names(5)=(/'In','Im','Is','C ','D '/) !names of situations
    integer, parameter :: inps(5)=(/0,1,0,0,0/)                 !initial number of people in each situation
    real,    parameter :: In(5)=(/0.5,0.5,0.,0.,0./)            !stochastic vector (In=Incubation)
    real,    parameter :: Im(5)=(/0.,0.6,0.3,0.1,0./)           !stochastic vector (Im=Infected mild)
    real,    parameter :: Is(5)=(/0.,0.,0.75,0.2,0.05/)         !stochastic vector (Is=Infected severe)
    real,    parameter :: C(5)=(/0.,0.,0.,1.,0./)               !stochastic vector (C=Cured)
    real,    parameter :: D(5)=(/0.,0.,0.,0.,1./)               !stochastic vector (D=Dead)

    !do not change these definitions
    integer, parameter   :: tc=floor(days*dc)       !total checks
    character(len=2)     :: space(x,y)              !space of individuals and their condition
    character(len=2)     :: cp(2)                   !condition of the person 1,2
    integer              :: pp(2,2)                 !position of person 1,2 (1=x,2=y)
    integer              :: alive                   !alive persons
    integer              :: contacts                !contacts per check
    integer              :: ci,cj,ck,cq,cw          !counters
    integer              :: random_integer          !function
    real                 :: ra_num                  !random number
    real                 :: corrector               !rounding corrector
    real                 :: SM(5,5)                 !Stochastic Matrix

    call init_random_seed

    open(1,file="Ovito.dat")

    !creating a stochastic matrix
    SM(1,:)=In
    SM(2,:)=Im
    SM(3,:)=Is
    SM(4,:)=C
    SM(5,:)=D

    !transformation of probabilities
    do cq=1,size(names)
        if (SM(cq,cq)==1.) go to 1
        SM(cq,cq)=dc*SM(cq,cq)/(1-SM(cq,cq)+dc*SM(cq,cq))
        ra_num=SM(cq,cq)
        sm(cq,:)=(1-SM(cq,cq))*sm(cq,:)/(sum(SM(cq,:))-SM(cq,cq))
        SM(cq,cq)=ra_num
1   end do

    !Monte Carlo calculations
    do cq=1,NMC
        print*, cq
        space='H'
        alive=x*y-inps(5)
        corrector=0.

        !random distribution of situations in the population
        do ci=1,size(names)
            do ck=1,inps(ci)
2               pp(1,:)=(/random_integer(1,x),random_integer(1,y)/)

                if (space(pp(1,1),pp(1,2))/='H') go to 2

                space(pp(1,1),pp(1,2))=names(ci)
            end do
        end do

        !Ovito visualization
        if ((cq==1).and.(visualization=='Y')) call Ovito(x,y,space,len(space),1)

        !running time
        do ck=1,tc
            !contacts number calculation
            contacts=nint(dcp*alive/(2*dc)+corrector)
            corrector=corrector+dcp*alive/(2*dc)-contacts

            !random contacts and infections
            do cj=1,contacts
3               call random_contact(x,y,distance,Lp,pp,cp,space,PBCs)

                !ban meeting with dead
                if ((cp(1)=='D').or.(cp(2)=='D')) goto 3

                !healthy-infected contact
                if ((cp(1)=='H').and.((cp(2)=='Im').or.(cp(2)=='Is'))) then
                    call random_number(ra_num)

                    if (ra_num<infectiveness) then
                        space(pp(1,1),pp(1,2))='In'
                    endif

                    elseif (((cp(1)=='Im').or.(cp(1)=='Is')).and.(cp(2)=='H')) then
                        call random_number(ra_num)

                        if (ra_num<infectiveness) then
                            space(pp(2,1),pp(2,2))='In'
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

                                    if (names(cw)=='D') alive=alive-1
                                    exit
                                endif
                            end do

                        case ('Im')
                            do cw=1,size(names)
                                if (ra_num<Sum(SM(2,1:cw))) then
                                    space(ci,cj)=names(cw)

                                    if (names(cw)=='D') alive=alive-1
                                    exit
                                endif
                            end do

                        case ('Is')
                            do cw=1,size(names)
                                if (ra_num<Sum(SM(3,1:cw))) then
                                    space(ci,cj)=names(cw)

                                    if (names(cw)=='D') alive=alive-1
                                    exit
                                endif
                            end do

                        case ('C')
                            do cw=1,size(names)
                                if (ra_num<Sum(SM(4,1:cw))) then
                                    space(ci,cj)=names(cw)

                                    if (names(cw)=='D') alive=alive-1
                                    exit
                                endif
                            end do

                        case ('D')
                            do cw=1,size(names)
                                if (ra_num<Sum(SM(5,1:cw))) then
                                    space(ci,cj)=names(cw)

                                    if (names(cw)=='D') alive=alive+1
                                    exit
                                endif
                            end do
                    end select
                end do
            end do

            if (alive<=1) exit

            !Ovito visualization
            if ((cq==1).and.(visualization=='Y').and.(mod(ck,tc/of)==0)) call Ovito(x,y,space,len(space),1)
        end do
    end do

    print*, 'END',nmc
end program Disease_spread_simulation


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
    integer,          intent(in)  :: x,y            !dimension lengths
    integer,          intent(in)  :: distance       !movement distance
    integer,          intent(out) :: pp(2,2)        !position of person 1,2 (1=x,2=y)
    character(len=2), intent(out) :: cp(2)          !condition of the person 1,2
    character(len=2), intent(in)  :: space(x,y)     !space of individuals and their condition
    character,        intent(in)  :: Lp             !shape of interaction, R=Rhombus, C=Circle, S=Square (Lp spaces)
    character,        intent(in)  :: PBCs           !periodic boundary conditions
    integer                       :: random_integer !function
    real                          :: n              !variable

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
end subroutine Ovito


subroutine init_random_seed()
    !random_seed based on system time
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
