program poisson_2d
    implicit none
  
    !----------------------------------------------------------------------
    ! Parâmetros do problema
    !----------------------------------------------------------------------
    integer, parameter :: Nx = 513        ! Número de pontos em x
    integer, parameter :: Ny = 513       ! Número de pontos em y
    real(8), parameter :: tol = 1.0e-6  ! Tolerância para o resíduo
  
    !----------------------------------------------------------------------
    ! Declaração de variáveis
    !----------------------------------------------------------------------
    real(8) :: dx, dy, res, resLocal
    real(8), dimension(Nx) :: x
    real(8), dimension(Ny) :: y
    real(8), dimension(Nx,Ny) :: u, unew, f
    integer :: i, j, iter, iCenter
    real(8) :: uExato
  
    !----------------------------------------------------------------------
    ! Definição do passo de malha e preenchimento dos vetores x e y
    !----------------------------------------------------------------------
    dx = 1.0d0 / (Nx - 1)
    dy = 1.0d0 / (Ny - 1)
  
    do i = 1, Nx
       x(i) = (i - 1) * dx
    end do
    do j = 1, Ny
       y(j) = (j - 1) * dy
    end do
  
    !----------------------------------------------------------------------
    ! Definição do termo-f (lado direito): f(x,y) = 2*(y - x)
    !----------------------------------------------------------------------
    do j = 1, Ny
       do i = 1, Nx
          f(i,j) = 2.d0 * ( y(j) - x(i) )
       end do
    end do
  
    !----------------------------------------------------------------------
    ! Inicialização da solução u e imposição da condição de contorno exata
    ! u_exato = y*(x^2 - x*y)
    !----------------------------------------------------------------------
    do j = 1, Ny
       do i = 1, Nx
          if ( i == 1 .or. i == Nx .or. j == 1 .or. j == Ny ) then
             u(i,j) = y(j)*( x(i)**2 - x(i)*y(j) )
          else
             u(i,j) = 0.d0  ! Chute inicial no interior
          end if
       end do
    end do
  
    unew = u  ! Cópia inicial
  
    !----------------------------------------------------------------------
    ! Abertura de arquivo para (iter, residuo)
    !----------------------------------------------------------------------
    open(unit=10, file="residuo.dat", status="unknown")
  
    !----------------------------------------------------------------------
    ! Iterações de Gauss-Seidel
    !    (u_{i+1,j} - 2u_{i,j} + u_{i-1,j}) / dx^2
    !  + (u_{i,j+1} - 2u_{i,j} + u_{i,j-1}) / dy^2 = f_{i,j}
    !
    ! ==> u_{i,j}^{(new)} =
    !           [ (u_{i+1,j}+u_{i-1,j})/dx^2 + (u_{i,j+1}+u_{i,j-1})/dy^2 - f(i,j) ]
    !           ---------------------------------------------------------------
    !                     [ 2/dx^2 + 2/dy^2 ]
    !----------------------------------------------------------------------
    res = 1.d10
    iter = 0
  
    do while ((res > tol))
       iter = iter + 1
       res = 0.d0
  
       ! Varre nós internos
       do j = 2, Ny - 1
          do i = 2, Nx - 1
             unew(i,j) = ( ( u(i+1,j) + u(i-1,j) ) / dx**2 + &
                           ( u(i,j+1) + u(i,j-1) ) / dy**2 -  &
                             f(i,j) )  &
                         / ( 2.d0 / dx**2 + 2.d0 / dy**2 )
          end do
       end do
  
       ! Calcula resíduo (norma L2 da diferença unew - u)
       do j = 2, Ny - 1
          do i = 2, Nx - 1
             resLocal = unew(i,j) - u(i,j)
             res = res + resLocal * resLocal
          end do
       end do
       res = sqrt( res / dble( (Nx-2)*(Ny-2) ) )
  
       ! Atualiza u
       u = unew
  
       ! Escreve (iter, res) no arquivo
       write(10,*) iter, res
    end do

   ! ...
   open(unit=30, file="solucao2d.dat", status="unknown")
   do j = 1, Ny
      do i = 1, Nx
         ! Grava X(i), Y(j) e U(i,j) em três colunas
         write(30,'(3E20.10)') x(i), y(j), u(i,j)
      end do
   end do
   close(30)

   print *, "Arquivo solucao2d.dat gerado com sucesso."
   ! ...

  
    close(10)
  
    !----------------------------------------------------------------------
    ! Escrita do corte vertical: x ~ 0.5  (iCenter corresponde a x=0.5)
    !----------------------------------------------------------------------
    iCenter = (Nx - 1) / 2 + 1
    open(unit=20, file="corte_vertical.dat", status="unknown")
  
    do j = 1, Ny
       uExato = y(j) * ( x(iCenter)**2 - x(iCenter)*y(j) )
       write(20,'(3E20.10)') y(j), u(iCenter,j), uExato
    end do
  
    close(20)
  
    !----------------------------------------------------------------------
    ! Mensagens finais
    !----------------------------------------------------------------------
    print *, "========================================"
    print *, "Concluída a iteração de Gauss-Seidel."
    print *, "Número de iterações =", iter
    print *, "Resíduo final       =", res
    print *, "Dados escritos em:"
    print *, " - residuo.dat (resíduo por iteração)"
    print *, " - corte_vertical.dat (solução numérica vs. exata em x=0.5)"
    print *, "========================================"
  
  end program poisson_2d
  