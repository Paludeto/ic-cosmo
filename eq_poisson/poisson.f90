! Condições de contorno de Dirchlet:
! L = W = 1 
! u(x, 0) = u(x, 1) = u(0, y) = u(1, y) = 0
! Equação discretizada para malhas uniformes: 
! u(i, j) = (u(i+1, j) + u(i-1, j) + u(i, j+1) + u(i, j-1)) / 4.0 - 2.0 * (y - x)


program poisson

    ! Declaração de variáveis
    implicit none
    integer :: i, j, iteracoes
    integer, parameter :: nX = 513, nY = 513
    real, parameter :: Lx = 1.0, Ly = 1.0, tol = 1.0e-6
    real :: dx, dy, erro, u_novo

    ! Matrizes (matriz sol. u, f(x, y))
    real, dimension(nX, nY) :: u, f

    ! Grid de pontos x, grid de pontos y
    REAL, DIMENSION(nX) :: x
    REAL, DIMENSION(nY) :: y

    ! Cálculo (dx = dy por ser uma malha numérica uniforme)
    dx = Lx / (nX - 1)
    dy = Ly / (nY - 1)

    ! Gera pontos da malha numérica a distância de 1 / 512 cada
    do i = 1, nX
        x(i) = (i - 1) * dx
    end do 
    do j = 1, nY
        y(j) = (j - 1) * dy
    end do 

    ! Inicializa a matriz u com 0 e calcula o termo fonte f(x, y) para cada ponto da malha
    u = 0.0
    do i = 1, nX
        do j = 1, nY
            f(i, j) = 2.0 * (y(j) - x(i))
        end do
    end do 
    
    erro = 1.0
    iteracoes = 0

    ! Gauss-Seidel
    do while (erro > tol)

        erro = 0.0

        do j = 2, nY - 1
            do i = 2, nX - 1
                u_novo = 0.25 * (u(i - 1, j) + u(i + 1, j) + u(i, j - 1) + u(i, j + 1) - dx**2 * f(i, j))
                erro = max(erro, abs(u_novo - u(i, j)))
                u(i, j) = u_novo
            end do 
        end do 

        iteracoes = iteracoes + 1

    end do 

    write(*, '(A, I6, A, ES15.8)') "Convergiu em ", iteracoes, " iterações, com erro de ", erro
    print *

    ! Geração do arquivo "matriz.dat"
    open(unit=10, file="matriz.dat", status="replace")
        do i = nY, 1, -1  ! Escrevendo de cima para baixo, convenção MatLab
            do j = 1, nX
                write(10, '(I4, I4, F10.5)') j, nY - i + 1, u(i, j)
            end do
            write(10, *)
        end do
    close(10)

end program poisson