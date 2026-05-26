program tridia

    implicit none

    integer :: n
    integer :: i
    integer :: unit_salida

    real(8), allocatable :: matriz_A(:,:)

    print *, "Ingresa n:"
    read(*,*) n

    if (n <= 0) then
        print *, "n debe ser positivo."
        stop
    end if

    allocate(matriz_A(n, n))

    matriz_A = 0.0

    do i = 1, n
        matriz_A(i, i) = 2.0
    end do

    do i = 1, n - 1
        matriz_A(i, i + 1) = -1.0
        matriz_A(i + 1, i) = -1.0
    end do

    unit_salida = 10
    open(unit=unit_salida, file="tridia.dat", status="replace", action="write")

    write(unit_salida, *) n

    do i = 1, n
        write(unit_salida, *) matriz_A(i, 1:n)
    end do

    close(unit_salida)

    print *, "Archivo tridia.dat generado correctamente."

end program tridia