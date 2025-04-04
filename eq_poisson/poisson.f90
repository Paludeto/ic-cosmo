! Condições de contorno de Dirchlet:
! L = W = 1 
! u(1, y) =  y - y², u(x, 1) =  x² - x

program poisson

    ! Declaração de variáveis
    implicit none
    integer :: i, j, iteracoes
    integer, parameter :: nX = 513, nY = 513
    real, parameter :: Lx = 1.0, Ly = 1.0, tol = 1.0e-6
    real(8) :: dx, dy, erro, u_novo

    ! Matrizes (matriz sol. u, f(x, y))
    real(8), dimension(nX, nY) :: u, f

    ! Grid de pontos x, grid de pontos y
    real(8), dimension(nX) :: x
    real(8), dimension(nY) :: y

    ! Cálculo (dx = dy por ser uma malha numérica uniforme)
    dx = Lx / (nX - 1)
    dy = Ly / (nY - 1)

    ! Gera pontos da malha numérica a distância de 1 / 512 cada
    do i = 1, nX
        x(i) = real(i - 1) * dx
    end do 
    do j = 1, nY
        y(j) = real(j - 1) * dy
    end do 

    ! Inicializa a matriz u com 0, calcula o termo fonte f(x, y) para cada ponto da malha 
    u = 0.0
    do i = 1, nX
        do j = 1, nY
            f(i, j) =  2.d0 * (y(j) - x(i))
        end do
    end do

    ! Imposição das condições de contorno de Dirichlet
    do j = 1, nY
        u(nX, j) = y(j) - y(j)**2 ! = y - y²
    end do
    do i = 1, nX
        u(i, nY) = x(i)**2 - x(i) ! = x² - x
    end do 

    erro = 1.0
    iteracoes = 0

    ! Cria o arquivo para armazenar o resíduo por iteração
    open(unit=10, file="matlab/dat/residuo.dat", status="replace")

    ! Gauss-Seidel
    do while (erro > tol)

        erro = 0.0

        do j = 2, nY - 1
            do i = 2, nX - 1
                u_novo = 0.25 * ( u(i+1,j) + u(i-1,j) + u(i,j+1) + u(i,j-1) - (dx**2) * f(i,j))
                erro = max(erro, abs(u_novo - u(i, j)))
                u(i, j) = u_novo
            end do 
        end do 

        iteracoes = iteracoes + 1

        ! Escreve iteração + erro atual em um arquivo .dat
        write(10, '(I6, ES25.8)') iteracoes, erro

    end do 

    close(10)

    write(*, '(A, I6, A, ES15.8)') "Convergiu em ", iteracoes, " iterações, com erro de ", erro
    print *

    ! Geração do arquivo "matriz.dat"
    open(unit=20, file="matlab/dat/matriz.dat", status="replace")
        do i = nY, 1, -1  ! Escrevendo de cima para baixo, convenção MatLab
            do j = 1, nX
                write(20, '(I4, I4, ES25.10E3)') j, nY - i + 1, u(i, j)
            end do
            write(20, *)
        end do
    close(20)

    ! Escrevendo a linha central (nX / 2) da matriz u em um .dat
    open(unit=30, file="matlab/dat/l_central.dat", status="replace")
        do j = nY, 1, -1
            write(30, '(F10.5, ES25.10E3)') real(j - 1) / real(nY - 1), u(int((nX / 2.0)) + 1, j)
        end do 
        write(30, *)
    close(30)

end program poisson