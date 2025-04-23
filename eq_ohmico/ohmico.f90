! Módulo que contém as constantes do problema
module constantes
    integer, parameter :: nX = 10, nY = 10, nZ = 10, max_iter = 10000
    real(8), parameter :: tol = 1.0d-5, L = 0.1d0, delta = L / (nX - 1)
    real(8), parameter :: B = 0.0118d0, C = 0.4888d0
    real(8), parameter :: tTotal = 1000.d0, nT = 1000.d0, dt = tTotal / nT
    real(8), parameter :: k_mat = 0.51d0, h_conv = 4.4d0, rho = 1070.0d0, cP = 3411.0d0, eta = 1.d0 / (rho * cP)
end module constantes

program ohmico
    use constantes
    implicit none

    real(8) :: T(nX, nY, nZ), V(nX, nY, nZ), U(nX, nY, nZ), T_new, erro, tempo_atual
    real(8) :: sigma_x_pos, sigma_x_ant, sigma_y_pos, sigma_y_ant, sigma_z_pos, sigma_z_ant
    integer :: i, j, k, iter

    ! Inicialização das matrizes
    T = 0.d0
    V = 0.d0
    U = 0.d0

    ! Condições de contorno
    do j = 1, nY
        do k = 1, nZ
            V(1,j,k) = 35.d0
            T(1,j,k) = 25.d0
        end do
    end do

    do i = 1, nX
        do k = 1, nZ
            T(i,1,k) = 25.d0
        end do
    end do

    tempo_atual = 0.d0

    do while (tempo_atual <= tTotal)

        erro = 1.d0
        iter = 0

        do while (erro > tol .and. iter < max_iter)
            erro = 0.d0
            iter = iter + 1

            do i = 2, nX - 1
                do j = 2, nY - 1
                    do k = 2, nZ - 1
                        sigma_x_pos = 0.5d0 * (condutividade(T(i,j,k)) + condutividade(T(i+1,j,k)))
                        sigma_x_ant = 0.5d0 * (condutividade(T(i,j,k)) + condutividade(T(i-1,j,k)))
                        sigma_y_pos = 0.5d0 * (condutividade(T(i,j,k)) + condutividade(T(i,j+1,k)))
                        sigma_y_ant = 0.5d0 * (condutividade(T(i,j,k)) + condutividade(T(i,j-1,k)))
                        sigma_z_pos = 0.5d0 * (condutividade(T(i,j,k)) + condutividade(T(i,j,k+1)))
                        sigma_z_ant = 0.5d0 * (condutividade(T(i,j,k)) + condutividade(T(i,j,k-1)))

                        V(i,j,k) = ( &
                            sigma_x_pos * V(i+1,j,k) + sigma_x_ant * V(i-1,j,k) + &
                            sigma_y_pos * V(i,j+1,k) + sigma_y_ant * V(i,j-1,k) + &
                            sigma_z_pos * V(i,j,k+1) + sigma_z_ant * V(i,j,k-1) ) / &
                            (sigma_x_pos + sigma_x_ant + sigma_y_pos + sigma_y_ant + sigma_z_pos + sigma_z_ant)
                    end do
                end do
            end do

            do i = 2, nX - 1
                do j = 2, nY - 1
                    do k = 2, nZ - 1
                        U(i,j,k) = eta * condutividade(T(i,j,k)) * ( &
                            ((V(i+1,j,k) - V(i-1,j,k)) / (2.d0 * delta)) ** 2 + &
                            ((V(i,j+1,k) - V(i,j-1,k)) / (2.d0 * delta)) ** 2 + &
                            ((V(i,j,k+1) - V(i,j,k-1)) / (2.d0 * delta)) ** 2 )
                    end do
                end do
            end do

            do i = 2, nX - 1
                do j = 2, nY - 1
                    do k = 2, nZ - 1
                        T_new = T(i,j,k) + (dt / (rho * cP)) * &
                            (k_mat * (T(i+1,j,k) + T(i-1,j,k) + T(i,j+1,k) + T(i,j-1,k) + T(i,j,k+1) + T(i,j,k-1) - 6.d0 * T(i,j,k)) / (delta ** 2) + U(i,j,k))

                        erro = max(erro, abs(T_new - T(i,j,k)))
                        T(i,j,k) = T_new
                    end do
                end do
            end do
        end do

        tempo_atual = tempo_atual + dt

    end do

    print *, "Algoritmo finalizado"

    open(unit=20, file="matlab/dat/T_volume.dat", status="replace")
    do k = 1, nZ
        do j = 1, nY
            do i = 1, nX
                write(20, *) i, j, k, T(i, j, k)
            end do
        end do
    end do
    close(20)

contains

    real(8) function condutividade(temperatura)
        real(8), intent(in) :: temperatura
        condutividade = B * temperatura + C
    end function

end program