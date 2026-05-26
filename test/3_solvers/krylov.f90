program krylov

    implicit none

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
    ! x_exacto
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

    integer :: n
    integer :: numero_maximo_krylov
    integer :: estado
    integer :: unidad_entrada

    real(8), allocatable :: matriz_A(:,:)
    real(8), allocatable :: vector_b(:)
    real(8), allocatable :: x_exacto(:)
    real(8), allocatable :: solucion_aproximada(:)
    real(8), allocatable :: residuo(:)

    real(8), allocatable :: matriz_krylov(:,:)
    real(8), allocatable :: matriz_krylov_actual(:,:)
    real(8), allocatable :: matriz_A_por_krylov(:,:)

    real(8), allocatable :: matriz_normal(:,:)
    real(8), allocatable :: vector_normal(:)
    real(8), allocatable :: coeficientes_krylov(:)

    real(8), allocatable :: vector_temporal(:)
    real(8) :: norma_vector_temporal
    real(8) :: norma_residuo
    real(8) :: norma_error
    real(8) :: factor
    real(8) :: suma

    integer :: i
    integer :: j
    integer :: k
    integer :: dimension_krylov

    unidad_entrada = 20
    open(unit=unidad_entrada, file="tridia.dat", status="old", action="read", iostat=estado)

    read(unidad_entrada, *, iostat=estado) n

    allocate(matriz_A(n, n))
    allocate(vector_b(n))
    allocate(x_exacto(n))
    allocate(solucion_aproximada(n))
    allocate(residuo(n))
    allocate(matriz_krylov(n, n))
    allocate(matriz_krylov_actual(n, n))
    allocate(matriz_A_por_krylov(n, n))
    allocate(matriz_normal(n, n))
    allocate(vector_normal(n))
    allocate(coeficientes_krylov(n))
    allocate(vector_temporal(n))

    do i = 1, n
        read(unidad_entrada, *, iostat=estado) matriz_A(i, 1:n)
        if (estado /= 0) then
            print *, "No se pudo leer la matriz desde tridia.dat."
            close(unidad_entrada)
            stop
        end if
    end do

    close(unidad_entrada)

    numero_maximo_krylov = n

    x_exacto = 1.0

    vector_b = matmul(matriz_A, x_exacto)

    matriz_krylov = 0.0
    matriz_krylov(:, 1) = vector_b

    print *, "Tamano del sistema:", n
    print *, " "

    do dimension_krylov = 1, numero_maximo_krylov

        if (dimension_krylov > 1) then

            vector_temporal = matmul(matriz_A, matriz_krylov(:, dimension_krylov - 1))

            norma_vector_temporal = sqrt(dot_product(vector_temporal, vector_temporal))

            matriz_krylov(:, dimension_krylov) = vector_temporal / norma_vector_temporal

        end if

        matriz_krylov_actual = 0.0
        matriz_A_por_krylov = 0.0
        matriz_normal = 0.0
        vector_normal = 0.0
        coeficientes_krylov = 0.0

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
!! Solver Gauss sin pivote de la normal
        do k = 1, dimension_krylov - 1
            do i = k + 1, dimension_krylov
                factor = matriz_normal(i, k) / matriz_normal(k, k)

                do j = k, dimension_krylov
                    matriz_normal(i, j) = matriz_normal(i, j) - factor * matriz_normal(k, j)
                end do

                vector_normal(i) = vector_normal(i) - factor * vector_normal(k)
            end do
        end do

        coeficientes_krylov(1:dimension_krylov) = 0.0

        do i = dimension_krylov, 1, -1
            suma = 0.0

            do j = i + 1, dimension_krylov
                suma = suma + matriz_normal(i, j) * coeficientes_krylov(j)
            end do

            coeficientes_krylov(i) = (vector_normal(i) - suma) / matriz_normal(i, i)
        end do
!! Fin solver
        solucion_aproximada = &
            matmul( &
                matriz_krylov_actual(:, 1:dimension_krylov), &
                coeficientes_krylov(1:dimension_krylov) &
            )

        residuo = vector_b - matmul(matriz_A, solucion_aproximada)

        norma_residuo = sqrt(dot_product(residuo, residuo))

        print '(A, I3, A, ES12.4)', "m = ", dimension_krylov, &
            "    norma residuo = ", norma_residuo

    end do

    norma_error = sqrt(dot_product(x_exacto - solucion_aproximada, &
                                  x_exacto - solucion_aproximada))

    print *, " "
    print *, "Xm Final:"
    do i = 1, n
        print '(I3, F14.8)', i, solucion_aproximada(i)
    end do

    print *, " "
    print *, "Error:"
    print *, norma_error

end program krylov