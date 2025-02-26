! Considerando uma placa retangular de dimensões L×W onde L = 1.0m e W = 1.0m.
! A placa está em regime estacionário de condução de calor, sem geração interna de calor. 
! Dessa forma, a distribuição de temperatura T(x,y) no interior da placa obedece à equação de Laplace:

! Condições de Contorno:
! As condições de contorno são fixadas da seguinte forma:
! Borda Superior: T(x,W) = 100∘C para 0 ≤ x ≤ L.
! Borda Inferior: T(x, 0) = 0∘C para 0 ≤ x ≤ L.
! Bordas Laterais: T(0,y) = 50∘C para 0 ≤ y ≤ W.

program laplace

    implicit none
    integer :: i, j, iteracoes, maxIter
    integer, parameter :: nX = 100, nY = 100
    real :: err, tol, maxErro, temp
    real :: L, W
    real, dimension(nX, nY) :: matriz = 0.0

    ! Parâmetros do método iterativo
    maxIter = 1000000000
    tol = 1.0e-12
    iteracoes = 0

    ! Malha uniforme
    L = 1.0
    W = 1.0 

    ! Borda superior
    do i = 1, nX
        matriz(i, nY) = 100.0
    end do

    ! Bordas laterais
    do j = 1, nY
        matriz(1, j) = 50.0
        matriz(nX, j) = 50.0
    end do

    ! Solução da equação de Laplace
    do
        maxErro = 0.0
        
        do j = 2, nY - 1
            do i = 2, nX - 1

                temp = matriz(i, j)

                matriz(i, j) = 0.25 * (matriz(i - 1, j) + matriz(i + 1, j) + matriz(i, j - 1) + matriz(i, j + 1))
                err = abs(matriz(i, j) - temp)

                if (err > maxErro) maxErro = err

            end do
        end do

        iteracoes = iteracoes + 1

        if (maxErro < tol) then
            print *, "Convergência atingida após", iteracoes, "iterações"
            exit
        end if

        if (iteracoes >= maxIter) then
            print *, "Não foi possível encontrar a solução em", iteracoes, "iterações"
            exit
        end if

    end do

    ! Impressão da matriz
    print *, "Matriz:"
    do j = nY, 1, -1
         write(*,*) (matriz(i, j), i = 1, nX)
    end do

    ! Geração do arquivo "matriz.dat"
    open(unit=10, file="matriz.dat", status="replace")
        do i = nY, 1, -1  ! Escrevendo de cima para baixo, convenção MatLab
            do j = 1, nX
                write(10, '(I4, I4, F10.5)') j, nY - i + 1, matriz(i, j)
            end do
            write(10, *)
        end do
    close(10)

end program laplace