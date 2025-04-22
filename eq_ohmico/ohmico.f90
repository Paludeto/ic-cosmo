! Módulo que contém as constantes do problema
module constantes
    integer, parameter:: nX = 10, nY = 10, nZ = 10, nT = 1000, tTotal = 1000, max_iter = 10000
    real(8), parameter :: tol = 1.0d-5, L = 0.1d0, delta = L / (nX - 1), dt = tTotal / nT
    real(8), parameter :: B = 0.0118d0, C = 0.4888d0
    real(8), parameter :: k_mat = 0.51d0, h_conv = 4.4d0, rho = 1070.0d0, cP = 3411.0d0, eta = 1.d0 / (rho * cP)
end module constantes

program ohmico

    use constantes
    implicit none

    ! Matriz para valores temperatura, potencial e energia
    real(8) :: T(nX, nY, nZ), V(nX, nY, nZ), U(nX, nY, nZ), T_new, erro
    real(8) :: sigma_x_pos, sigma_x_ant, sigma_y_pos, sigma_y_ant, sigma_z_pos, sigma_z_ant
    integer :: i, j, k, iter

    ! Preenchendo matrizes com zeros (T_inf no caso da temperatura)
    T = 0.d0
    V = 0.d0
    U = 0.d0

    ! Condições de contorno
    ! Plano YZ
    do j = 1, nY
        do k = 1, nZ
            V(1,j,k) = 35.d0        
            T(1,j,k) = 25.d0
        end do
    end do

    ! Plano XZ
    do i = 1, nX
        do k = 1, nZ
            T(i,1,k) = 25.d0         
        end do
    end do

    erro = 1.d0
    iter = 0

    ! Adicionar transiência, variação de acordo com o tempo, plotar queda resíduo, adicionar cond contorno em V e T.
    do while (erro > tol .and. iter < max_iter) 

        erro = 0.d0
        iter = iter + 1

        do i = 2, nX - 1
            do j = 2, nY - 1
                do k = 2, nZ - 1

                    ! Sigma X
                    sigma_x_pos = (condutividade(T(i, j, k)) + condutividade(T(i + 1, j, k))) * 0.5d0   
                    sigma_x_ant = (condutividade(T(i, j, k)) + condutividade(T(i - 1, j, k))) * 0.5d0

                    ! Sigma Y
                    sigma_y_pos = (condutividade(T(i, j, k)) + condutividade(T(i, j + 1, k))) * 0.5d0
                    sigma_y_ant = (condutividade(T(i, j, k)) + condutividade(T(i, j - 1, k))) * 0.5d0

                    ! Sigma Z
                    sigma_z_pos = (condutividade(T(i, j, k)) + condutividade(T(i, j, k + 1))) * 0.5d0
                    sigma_z_ant = (condutividade(T(i, j, k)) + condutividade(T(i, j, k - 1))) * 0.5d0

                    ! Matriz potencial elétrico
                    V(i,j,k) = ( &
                        sigma_x_pos * V(i + 1, j, k) + sigma_x_ant * V(i - 1, j, k) + &
                        sigma_y_pos * V(i, j + 1, k) + sigma_y_ant * V(i, j - 1,k) + &
                        sigma_z_pos * V(i, j, k + 1) + sigma_z_ant * V(i, j, k - 1) ) / &
                        (sigma_x_pos + sigma_x_ant + sigma_y_pos + sigma_y_ant + sigma_z_pos + sigma_z_ant)
                                
                end do 
            end do 
        end do 

        do i = 2, nX - 1
            do j = 2, nY - 1
                do k = 2, nZ - 1

                    U(i, j, k) = eta * condutividade(T(i, j, k)) * ( &
                    ((V(i + 1, j, k) - V(i - 1, j, k)) / (2.d0 * delta)) ** 2 + &
                    ((V(i, j + 1, k) - V(i, j - 1, k)) / (2.d0 * delta)) ** 2 + &
                    ((V(i, j, k + 1) - V(i, j, k - 1)) / (2.d0 * delta)) ** 2 )
                    
                end do
            end do 
        end do 
    
        do i = 2, nX - 1
            do j = 2, nY - 1
                do k = 2, nZ - 1
                    
                    T_new = T(i, j, k) + (dt / (rho * cP)) * &
                    (k_mat * (T(i + 1, j, k) + T(i - 1, j, k) + &
                        T(i, j + 1, k) + T(i, j - 1, k) + &
                        T(i, j, k+1) + T(i, j, k - 1) - 6.d0 * T(i, j, k) &
                    ) / (delta ** 2) + U(i, j, k))
    
                    erro = max(erro, abs(T_new - T(i,j,k)))
                    T(i, j, k) = T_new

                end do 
            end do 
        end do 

    end do 

    print *,  "Algoritmo finalizado em", iter, "iterações"

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

    ! 𝜎(𝑇)= 𝐵𝑇 + C
    ! Função de condutividade, armazena resultado no parâmetro resultado
    real(8) function condutividade(temperatura)

        real(8), intent(in) :: temperatura

        condutividade = (B * temperatura) + C

    end function


end program