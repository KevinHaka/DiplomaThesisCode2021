Program SIR_DSS
    implicit none

    !initial conditions and parameters
    integer,   parameter :: NMC=10           !Number Monte Carlo
    integer,   parameter :: x=100,y=100       !dimension lengths (x*y=population)
    integer,   parameter :: ii=50             !initial infected
    integer,   parameter :: ir=0              !initial removed
    integer,   parameter :: of=100            !Ovito frames (of<day*dc)
    real,      parameter :: days=50.          !time
    real,      parameter :: dc=50.            !daily checks (dc->infinity => dt->0)
    real,      parameter :: g=0.05            !removal rate (g<=dc)
    real,      parameter :: b=0.00004         !infection rate
    character, parameter :: visualization='N' !Ovito visualization results (Y=Yes,N=No)

    !do not change these definitions
    integer, parameter :: tc=floor(days*dc)   !total checks
    integer, parameter :: is=x*y-ii-ir        !initial susceptible
    real,    parameter :: P_I=1-g/dc          !Probability of staying in situation I
    character          :: space(x,y)          !space of individuals and their condition
    character          :: cp(2)               !condition of the person 1,2
    integer            :: pp(2,2)             !position of person 1,2 (1=x,2=y)
    integer            :: nSIR(3)             !number of people in each situation (1=S:Susceptible, 2=I:Infected, 3=R:Removed)
    integer            :: contacts            !contacts in population (healthy-infected)
    integer            :: ci,cj,ck,cq         !counters
    integer            :: random_integer      !function
    real               :: ra_num              !random number
    real               :: SIR_MCc(3,tc+1)     !SIR Monte Carlo calculations
    real               :: corrector           !rounding corrector

    call init_random_seed

    open(1,file="Ovito.dat")
    open(2,file="Results.dat")

    !Monte Carlo calculations
    SIR_MCc=0.
    do cq=1,NMC
        print*, cq
        space='S'
        corrector=0.
        nSIR=(/is,ii,ir/)

        !random infections in the population
        do ck=1,ii
1           pp(1,:)=(/random_integer(1,x),random_integer(1,y)/)

            if (space(pp(1,1),pp(1,2))=='S') then
                space(pp(1,1),pp(1,2))='I'

                else
                    goto 1
            end if
        end do

        !random removes in the population
        do ck=1,ir
2           pp(1,:)=(/random_integer(1,x),random_integer(1,y)/)

            if (space(pp(1,1),pp(1,2))=='S') then
                space(pp(1,1),pp(1,2))='R'

                else
                    goto 2
            end if
        end do

        SIR_MCc(:,1)=SIR_MCc(:,1)+nSIR/Real(NMC)

        !Ovito visualization
        if ((cq==1).and.(visualization=='Y')) call Ovito(x,y,space,1)

        !running time
        do ck=1,tc
            !random contacts and infections
            contacts=nint((b*nSIR(1)*nSIR(2))/dc+corrector)
            corrector=corrector+(b*nSIR(1)*nSIR(2))/dc-contacts

            do cj=1,contacts
                !safety condition
                if ((nSIR(1)==0).or.(nSIR(2)==0)) exit

3               pp(1,:)=(/random_integer(1,x),random_integer(1,y)/)
                pp(2,:)=(/random_integer(1,x),random_integer(1,y)/)
                cp=(/space(pp(1,1),pp(1,2)),space(pp(2,1),pp(2,2))/)

                if (((cp(1)=='S').and.(cp(2)=='I')).or.((cp(1)=='I').and.(cp(2)=='S'))) then
                    selectcase (cp(1))
                        case ('S')
                            space(pp(1,1),pp(1,2))='I'
                            nSIR(1)=nSIR(1)-1
                            nSIR(2)=nSIR(2)+1

                        case ("I")
                            space(pp(2,1),pp(2,2))='I'
                            nSIR(1)=nSIR(1)-1
                            nSIR(2)=nSIR(2)+1
                    endselect

                    else
                        goto 3
                end if
            end do

            SIR_MCc(:,ck+1)=SIR_MCc(:,ck+1)+nSIR/Real(NMC)

            !disease progression in each individual
            do ci=1,x
                do cj=1,y
                    if (space(ci,cj)=='I') then
                        call random_number(ra_num)

                        if (ra_num>P_I) then
                            space(ci,cj)='R'
                            nSIR(2)=nSIR(2)-1
                            nSIR(3)=nSIR(3)+1
                        end if
                    end if
                end do
            end do

            !Ovito visualization
            if ((cq==1).and.(visualization=='Y').and.(mod(ck,tc/of)==0)) call Ovito(x,y,space,1)
        end do
    enddo

    !writes the results
    do ci=1,tc+1
        if (mod(ci,int(dc))==1) write(2,*) (ci-1)/dc,SIR_MCc(:,ci)
    end do
end program SIR_DSS


integer function random_integer(a,b)
    !random integer value
    implicit none
    real                :: ra_var !random variable
    integer, intent(in) :: a,b    !space limits

    call random_number(ra_var)
    random_integer=floor((b-a+1)*ra_var)+a
end function


subroutine Ovito(x,y,space,nn)
    !Ovito visualization
    implicit none
    integer,          intent(in) :: nn         !file number (open in the main program)
    integer,          intent(in) :: x,y        !dimension lengths
    character(len=1), intent(in) :: space(x,y) !space of individuals and their condition
    integer                      :: i,j        !counters

    write(nn,*) x*y
    write(nn,*)
    do i=1,x
        do j=1,y
            write(nn,*) space(i,j),i,j
        end do
    end do
end subroutine


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
