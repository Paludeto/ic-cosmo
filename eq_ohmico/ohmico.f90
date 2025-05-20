module constantes
    implicit none
    integer, parameter :: nX=30, nY=30, nZ=30, max_iter=10000
    real(8), parameter :: tol=1.0d-5, L=0.1d0
    real(8), parameter :: deltaX=L/(nX-1), deltaY=L/(nY-1), deltaZ=L/(nZ-1)
    real(8), parameter :: B=0.0118d0, C=0.4888d0
    real(8), parameter :: tTotal=1000.d0, nT=1000.d0, dt=tTotal/nT
    real(8), parameter :: k_mat=0.51d0, h_conv=4.4d0
    real(8), parameter :: rho=1070.d0, cP=3411.d0, eta=1.d0/(rho*cP)
    real(8), parameter :: T_inf=35.d0      ! temperatura ambiente
    real(8), parameter :: V_applied=1.0d0    ! potencial aplicado (ajuste se precisar)
    integer, parameter :: j_start = nY/3 + 1, j_end = nY - nY/3
    integer, parameter :: k_start = nZ/3 + 1, k_end = nZ - nZ/3
end module constantes

program ohmico
    use constantes
    implicit none

    real(8) :: T(nX,nY,nZ), V(nX,nY,nZ), U(nX,nY,nZ)
    real(8) :: T_new, erro, tempo_atual
    real(8) :: sxp, sax, syp, say, szp, saz
    integer :: i,j,k,iter

    T = 0;  V = 0;  U = 0
    tempo_atual = 0

    do while (tempo_atual <= tTotal)
        erro = 1;  iter = 0
        do while (erro > tol .and. iter < max_iter)
            iter = iter + 1

            !––– Condição de contorno para V –––
            ! faces x=0 e x=L com square central (Dirichlet) e Neumann fora
            do j=1,nY; do k=1,nZ
                if (j>=j_start .and. j<=j_end .and. k>=k_start .and. k<=k_end) then
                    V(1 ,j,k) = V_applied
                    V(nX,j,k) = 0.d0
                else
                    V(1 ,j,k) = V(2 ,j,k)
                    V(nX,j,k) = V(nX-1,j,k)
                endif
            end do; end do
            ! Neumann nas demais faces
            do i=1,nX; do k=1,nZ
                V(i,1 ,k) = V(i,2 ,k)
                V(i,nY,k) = V(i,nY-1,k)
            end do; end do
            do i=1,nX; do j=1,nY
                V(i,j,1 ) = V(i,j,2 )
                V(i,j,nZ) = V(i,j,nZ-1)
            end do; end do

            !––– Gauss-Seidel interior para V –––
            do i=2,nX-1; do j=2,nY-1; do k=2,nZ-1
                sxp = 0.5d0*(condutividade(T(i,j,k))+condutividade(T(i+1,j,k)))
                sax = 0.5d0*(condutividade(T(i,j,k))+condutividade(T(i-1,j,k)))
                syp = 0.5d0*(condutividade(T(i,j,k))+condutividade(T(i,j+1,k)))
                say = 0.5d0*(condutividade(T(i,j,k))+condutividade(T(i,j-1,k)))
                szp = 0.5d0*(condutividade(T(i,j,k))+condutividade(T(i,j,k+1)))
                saz = 0.5d0*(condutividade(T(i,j,k))+condutividade(T(i,j,k-1)))
                V(i,j,k) = ( &
                  sxp*V(i+1,j,k)+sax*V(i-1,j,k)+ &
                  syp*V(i,j+1,k)+say*V(i,j-1,k)+ &
                  szp*V(i,j,k+1)+saz*V(i,j,k-1) ) / &
                  (sxp+sax+syp+say+szp+saz)
            end do; end do; end do

            !––– Atualiza U interior –––
            do i=2,nX-1; do j=2,nY-1; do k=2,nZ-1
                U(i,j,k) = eta*condutividade(T(i,j,k)) * &
                  (((V(i+1,j,k)-V(i-1,j,k))/(2.d0*deltaX))**2 + &
                   ((V(i,j+1,k)-V(i,j-1,k))/(2.d0*deltaY))**2 + &
                   ((V(i,j,k+1)-V(i,j,k-1))/(2.d0*deltaZ))**2)
            end do; end do; end do

            !––– Atualiza T interior –––
            erro = 0.d0
            do i=2,nX-1; do j=2,nY-1; do k=2,nZ-1
                T_new = T(i,j,k) + dt/(rho*cP) * &
                  ( k_mat*( (T(i+1,j,k)-2.d0*T(i,j,k)+T(i-1,j,k))/(deltaX**2) + &
                            (T(i,j+1,k)-2.d0*T(i,j,k)+T(i,j-1,k))/(deltaY**2) + &
                            (T(i,j,k+1)-2.d0*T(i,j,k)+T(i,j,k-1))/(deltaZ**2) ) + &
                    U(i,j,k) )
                erro = max(erro, abs(T_new - T(i,j,k)))
                T(i,j,k) = T_new
            end do; end do; end do

            !––– Condição de contorno para T –––
            do j=1,nY; do k=1,nZ
                if (j>=j_start .and. j<=j_end .and. k>=k_start .and. k<=k_end) then
                    ! Robin central
                    T(1 ,j,k) = ((k_mat/deltaX)*T(2 ,j,k) + h_conv*T_inf) / (k_mat/deltaX + h_conv)
                    T(nX,j,k) = ((k_mat/deltaX)*T(nX-1,j,k) + h_conv*T_inf) / (k_mat/deltaX + h_conv)
                else
                    ! Neumann fora
                    T(1 ,j,k) = T(2 ,j,k)
                    T(nX,j,k) = T(nX-1,j,k)
                endif
            end do; end do
            do i=1,nX; do k=1,nZ
                T(i,1 ,k) = T(i,2 ,k)
                T(i,nY,k) = T(i,nY-1,k)
            end do; end do
            do i=1,nX; do j=1,nY
                T(i,j,1 ) = T(i,j,2 )
                T(i,j,nZ) = T(i,j,nZ-1)
            end do; end do

        end do  ! fim GS interno

        tempo_atual = tempo_atual + dt
    end do      ! loop de tempo

    print *, "Algoritmo finalizado"
    open(unit=20,file="matlab/dat/T_volume.dat",status="replace")
      do k=1,nZ; do j=1,nY; do i=1,nX
        write(20,*) i,j,k,T(i,j,k)
      end do; end do; end do
    close(20)

contains

    real(8) function condutividade(temperatura)
        real(8), intent(in) :: temperatura
        condutividade = B*temperatura + C
    end function

end program