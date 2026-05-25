program krylov_simple

    implicit none

    ! ============================================================
    ! GLOSARIO
    !
    ! n
    !   Tamano del sistema Ax = b.
    !
    ! numero_maximo_krylov
    !   Dimension maxima del subespacio de Krylov.
    !
    ! matriz_A
    !   Matriz del sistema Ax = b.
    !
    ! vector_b
    !   Lado derecho del sistema.
    !
    ! solucion_exacta
    !   Solucion conocida artificialmente para verificar el codigo.
    !
    ! matriz_krylov
    !   Matriz cuyas columnas son b, Ab, A^2 b, ...
    !
    ! matriz_A_por_krylov
    !   Producto A * K_m.
    !
    ! matriz_normal
    !   Matriz pequena (A K_m)^T (A K_m).
    !
    ! vector_normal
    !   Vector pequeno (A K_m)^T b.
    !
    ! coeficientes_krylov
    !   Vector z que combina las columnas de Krylov.
    !
    ! solucion_aproximada
    !   Aproximacion x_m = K_m z.
    !
    ! residuo
    !   Vector r = b - A x_m.
    ! ============================================================

    integer, parameter :: n = 10
    integer, parameter :: numero_maximo_krylov = 10

    real(8) :: matriz_A(n, n)
    real(8) :: vector_b(n)
    real(8) :: solucion_exacta(n)
    real(8) :: solucion_aproximada(n)
    real(8) :: residuo(n)

    real(8) :: matriz_krylov(n, numero_maximo_krylov)
    real(8) :: matriz_krylov_actual(n, numero_maximo_krylov)
    real(8) :: matriz_A_por_krylov(n, numero_maximo_krylov)

    real(8) :: matriz_normal(numero_maximo_krylov, numero_maximo_krylov)
    real(8) :: vector_normal(numero_maximo_krylov)
    real(8) :: coeficientes_krylov(numero_maximo_krylov)

    real(8) :: vector_temporal(n)
    real(8) :: norma_vector_temporal
    real(8) :: norma_residuo
    real(8) :: norma_error

    integer :: i
    integer :: j
    integer :: dimension_krylov

    call construir_matriz_tridiagonal(n, matriz_A)

    solucion_exacta = 1.0d0

    vector_b = matmul(matriz_A, solucion_exacta)

    matriz_krylov = 0.0d0
    matriz_krylov(:, 1) = vector_b

    print *, "Metodo simple con subespacios de Krylov"
    print *, "Sistema Ax = b"
    print *, "Tamano del sistema:", n
    print *, " "

    do dimension_krylov = 1, numero_maximo_krylov

        if (dimension_krylov > 1) then

            vector_temporal = matmul(matriz_A, matriz_krylov(:, dimension_krylov - 1))

            norma_vector_temporal = sqrt(dot_product(vector_temporal, vector_temporal))

            matriz_krylov(:, dimension_krylov) = vector_temporal / norma_vector_temporal

        end if

        matriz_krylov_actual = 0.0d0
        matriz_A_por_krylov = 0.0d0
        matriz_normal = 0.0d0
        vector_normal = 0.0d0
        coeficientes_krylov = 0.0d0

        matriz_krylov_actual(:, 1:dimension_krylov) = matriz_krylov(:, 1:dimension_krylov)

        matriz_A_por_krylov(:, 1:dimension_krylov) = &
            matmul(matriz_A, matriz_krylov_actual(:, 1:dimension_krylov))

        matriz_normal(1:dimension_krylov, 1:dimension_krylov) = &
            matmul( &
                transpose(matriz_A_por_krylov(:, 1:dimension_krylov)), &
                matriz_A_por_krylov(:, 1:dimension_krylov) &
            )

        vector_normal(1:dimension_krylov) = &
            matmul( &
                transpose(matriz_A_por_krylov(:, 1:dimension_krylov)), &
                vector_b &
            )

        call resolver_sistema_pequeno( &
            dimension_krylov, &
            matriz_normal(1:dimension_krylov, 1:dimension_krylov), &
            vector_normal(1:dimension_krylov), &
            coeficientes_krylov(1:dimension_krylov) &
        )

        solucion_aproximada = &
            matmul( &
                matriz_krylov_actual(:, 1:dimension_krylov), &
                coeficientes_krylov(1:dimension_krylov) &
            )

        residuo = vector_b - matmul(matriz_A, solucion_aproximada)

        norma_residuo = sqrt(dot_product(residuo, residuo))

        print '(A, I3, A, ES12.4)', "dimension Krylov = ", dimension_krylov, &
            "    norma residuo = ", norma_residuo

    end do

    norma_error = sqrt(dot_product(solucion_exacta - solucion_aproximada, &
                                  solucion_exacta - solucion_aproximada))

    print *, " "
    print *, "Solucion aproximada final:"
    do i = 1, n
        print '(I3, F14.8)', i, solucion_aproximada(i)
    end do

    print *, " "
    print *, "Norma del error contra solucion exacta:"
    print *, norma_error

contains

    subroutine construir_matriz_tridiagonal(n, matriz_A)

        implicit none

        integer, intent(in) :: n
        real(8), intent(out) :: matriz_A(n, n)

        integer :: i

        matriz_A = 0.0d0

        do i = 1, n
            matriz_A(i, i) = 2.0d0
        end do

        do i = 1, n - 1
            matriz_A(i, i + 1) = -1.0d0
            matriz_A(i + 1, i) = -1.0d0
        end do

    end subroutine construir_matriz_tridiagonal


    subroutine resolver_sistema_pequeno(n, matriz_A, vector_b, vector_x)

        implicit none

        integer, intent(in) :: n
        real(8), intent(inout) :: matriz_A(n, n)
        real(8), intent(inout) :: vector_b(n)
        real(8), intent(out) :: vector_x(n)

        real(8) :: factor
        real(8) :: suma

        integer :: i
        integer :: j
        integer :: k

        do k = 1, n - 1
            do i = k + 1, n
                factor = matriz_A(i, k) / matriz_A(k, k)

                do j = k, n
                    matriz_A(i, j) = matriz_A(i, j) - factor * matriz_A(k, j)
                end do

                vector_b(i) = vector_b(i) - factor * vector_b(k)
            end do
        end do

        vector_x = 0.0d0

        do i = n, 1, -1
            suma = 0.0d0

            do j = i + 1, n
                suma = suma + matriz_A(i, j) * vector_x(j)
            end do

            vector_x(i) = (vector_b(i) - suma) / matriz_A(i, i)
        end do

    end subroutine resolver_sistema_pequeno

end program krylov_simple