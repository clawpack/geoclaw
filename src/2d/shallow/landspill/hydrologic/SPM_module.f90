!
! SPM_module.f90
! Copyright (C) 2018 Pi-Yueh Chuang <pychuang@gwu.edu>
!
! Distributed under terms of the MIT license.
!

!> @brief Naive sparse-matrix implementations.
module SPM_module
    implicit none
    private
    public:: IndexSet, COO, CSR, compress

    !> @brief A helper type representing indices and value of a non-zero for COO.
    type:: IndexSet
        private
        !> @brief Row and column index.
        integer(kind=4), dimension(2):: idx = [-999, -999]
        !> @brief Value.
        real(kind=8):: val = 0D0
        !> @brief Only used when being treated as linked list element in COO.
        type(IndexSet), pointer:: prev => null()
        !> @brief Only used when being treated as linked list element in COO.
        type(IndexSet), pointer:: next => null()

        contains
        !> @brief Initialization.
        procedure, private:: init_IndexSet
        !> @brief Output function.
        procedure, private:: write_IndexSet
        !> @brief Overloading of intrinsic write.
        generic:: write(formatted) => write_IndexSet
        !> @brief Destructor
        final:: destructor_IndexSet
    end type IndexSet

    !> @brief A mimic to C++ style constructor.
    interface IndexSet
        procedure:: constructor_IndexSet
    end interface IndexSet

    !> brief A naive implementation of COO sparse matrix
    type:: COO
        private
        !> @brief Number of non-zeros.
        integer(kind=4):: nnz = 0
        !> @brief Number of rows in the matrix.
        integer(kind=4):: n_rows = 0
        !> @brief Number of columns in the matrix.
        integer(kind=4):: n_cols = 0
        !> @brief Begining of the list.
        type(IndexSet), pointer:: bg => null()
        !> @brief The last element of the list.
        type(IndexSet), pointer:: ed => null()

        contains
        !> @brief Initialization.
        procedure:: init => init_COO
        !> @brief  Add a non-zero to the matrix.
        procedure:: append => append_COO
        !> @brief Destroy this COO.
        procedure:: destroy => destroy_COO
        !> @brief Output function.
        procedure, private:: write_COO
        !> @brief Overloading of intrinsic write.
        generic:: write(formatted) => write_COO
        !> @brief Destructor
        final:: destructor_COO
    end type COO

    !> @brief A mimic to C++ style constructor.
    interface COO
        procedure:: constructor_COO
    end interface COO

    !> @brief A naive implementation of CSR sparse matrix.
    type:: CSR
        private
        !> @brief Number of non-zeros.
        integer(kind=4):: nnz = 0
        !> @brief Number of rows in the matrix.
        integer(kind=4):: n_rows = 0
        !> @brief Number of columns in the matrix.
        integer(kind=4):: n_cols = 0
        !> @brief Number of non-zeros in each row.
        integer(kind=4), allocatable, dimension(:):: rows
        !> @brief Column indices.
        integer(kind=4), allocatable, dimension(:):: cols
        !> @brief Values of non-zeros.
        real(kind=8), allocatable, dimension(:):: vals

        contains
        !> @brief Init from a CSR matrix.
        procedure:: init => init_CSR
        !> @brief Access matrix element.
        procedure:: get => get_CSR
        !> @brief Set matrix element.
        procedure:: set => set_CSR
        !> @brief Add value to a matrix element.
        procedure:: add => add_CSR
        !> @brief Count nnz in a given block region.
        procedure:: count_block => count_block_CSR
        !> @brief Add values to non-zeros inside a given block region.
        procedure:: add_block => add_block_CSR
        !> @brief Destroy this CSR.
        procedure:: destroy => destroy_CSR
        !> @brief Write.
        procedure, private:: write_CSR
        !> @brief Overloading intrinsic write.
        generic:: write(formatted) => write_CSR
        !> @brief Destructor.
        final:: destructor_CSR
    end type CSR

    !> @brief An intergace to standard C qsort function
    interface
        subroutine qsort(arry, n_elem, size_of_elem, compare) bind(C, name='qsort')
            use, intrinsic:: iso_c_binding, only: c_size_t, c_int, c_ptr
            import
            type(*), intent(inout):: arry(*)
            integer(kind=c_size_t), intent(in), value:: n_elem
            integer(kind=c_size_t), intent(in), value:: size_of_elem
            abstract interface
                function compare_func(a, b) bind(C)
                    import
                    type(c_ptr), intent(in), value:: a, b
                    integer(kind=c_int):: compare_func
                end function compare_func
            end interface
            procedure(compare_func):: compare
        end subroutine qsort
    end interface

    !> @brief An intergace to standard C bsearch function
    interface
        function bsearch(key, arry, n_elem, size_of_elem, compare) bind(C, name='bsearch')
            use, intrinsic:: iso_c_binding, only: c_size_t, c_int, c_ptr
            import
            type(*), intent(in):: key
            type(*), intent(in):: arry(*)
            integer(kind=c_size_t), intent(in), value :: n_elem
            integer(kind=c_size_t), intent(in), value :: size_of_elem
            abstract interface
                function compare_func(a, b) bind(C)
                    import
                    type(c_ptr), intent(in), value :: a, b
                    integer(kind=c_int):: compare_func
                end function compare_func
            end interface
            procedure(compare_func):: compare
            type(c_ptr) :: bsearch
        end function bsearch
    end interface

contains

    ! constructor_IndexSet
    function constructor_IndexSet(row, col, val) result(IS)
        integer(kind=4), intent(in):: row, col
        real(kind=8), intent(in):: val
        type(IndexSet):: IS

        call init_IndexSet(IS, row, col, val)
    end function constructor_IndexSet
    
    ! init_IndexSet
    subroutine init_IndexSet(this, row, col, val)
        class(IndexSet), intent(inout):: this
        integer(kind=4), intent(in):: row, col
        real(kind=8), intent(in):: val

        this%idx = [row, col]
        this%val = val
    end subroutine init_IndexSet

    ! write_IndexSet
    subroutine write_IndexSet(this, iounit, iotype, v_list, stat, msg)
        ! variable declaration
        class(IndexSet), intent(in):: this
        integer(kind=4), intent(in):: iounit
        character(*), intent(in)::iotype
        integer(kind=4), intent(in):: v_list(:)
        integer(kind=4), intent(out):: stat
        character(*), intent(inout):: msg

        write(iounit, *, iostat=stat, iomsg=msg) &
            "(", this%idx(1), ", ", this%idx(2), "): ", this%val

    end subroutine write_IndexSet

    ! destructor_IndexSet
    subroutine destructor_IndexSet(this)
        type(IndexSet), intent(inout):: this
        this%idx = [-999, -999]
        this%val = 0D0
    end subroutine destructor_IndexSet

    ! constructor_COO
    function constructor_COO(n_rows, n_cols) result(mtx)
        integer(kind=4), intent(in):: n_rows, n_cols
        type(COO):: mtx

        call init_COO(mtx, n_rows, n_cols)
    end function constructor_COO
    
    ! init_COO
    subroutine init_COO(this, n_rows, n_cols)
        class(COO), intent(inout):: this
        integer(kind=4), intent(in):: n_rows, n_cols

        this%n_rows = n_rows
        this%n_cols = n_cols
    end subroutine init_COO

    ! append_COO
    subroutine append_COO(this, is)
        class(COO), intent(inout):: this
        type(IndexSet), intent(in):: is

        if (associated(this%ed)) then
            allocate(this%ed%next)
            this%ed%next%prev => this%ed
            this%ed => this%ed%next
        else ! the fisrt element in this COO
            allocate(this%bg)
            this%ed => this%bg
        end if

        this%ed%idx = is%idx
        this%ed%val = is%val
        this%nnz = this%nnz + 1
    end subroutine append_COO

    ! write_COO
    subroutine write_COO(this, iounit, iotype, v_list, stat, msg)
        ! variable declaration
        class(COO), intent(in):: this
        integer(kind=4), intent(in):: iounit
        character(*), intent(in)::iotype
        integer(kind=4), intent(in):: v_list(:)
        integer(kind=4), intent(out):: stat
        character(*), intent(inout):: msg
        integer(kind=4):: i
        character:: n, t
        type(IndexSet), pointer:: head


        n = new_line(t) ! n is the "new line" character
        t = achar(9) ! t is the character for a tab

        write(iounit, *, iostat=stat, iomsg=msg) "NNZ = ", this%nnz
        write(iounit, *, iostat=stat, iomsg=msg) n

        head => this%bg
        do i = 1, this%nnz
            write(iounit, *, iostat=stat, iomsg=msg) i, ":", &
                "(", head%idx(1), ",", head%idx(2), ") =>", head%val
            write(iounit, *, iostat=stat, iomsg=msg) n
            head => head%next
        end do
        nullify(head)
    end subroutine write_COO

    ! destroy_COO
    subroutine destroy_COO(this)
        class(COO), intent(inout):: this
        type(IndexSet), pointer:: head
        integer(kind=4):: i

        this%n_rows = 0
        this%n_cols = 0

        ! this means the list is not initialized
        if (.not. associated(this%bg)) return

        do i = this%nnz, 2, -1
            this%ed => this%ed%prev
            deallocate(this%ed%next)
        end do

        if (.not. associated(this%bg, this%ed)) then
            print *, "Error in deallocating COO."
            stop
        end if

        deallocate(this%bg)
        this%nnz = 0

    end subroutine destroy_COO

    ! destructor_COO
    subroutine destructor_COO(this)
        type(COO), intent(inout):: this
        call this%destroy()
    end subroutine destructor_COO

    ! constructor_CSR
    function constructor_CSR(src) result(mtx)
        type(COO), intent(in):: src
        type(CSR):: mtx

        call init_CSR(mtx, src)
    end function constructor_CSR

    ! init_CSR
    subroutine init_CSR(this, src)
        use, intrinsic:: iso_c_binding, only: c_size_t, c_sizeof
        class(CSR), intent(inout):: this
        type(COO), intent(in):: src
        integer(kind=4):: i, j, i_nnz
        type(IndexSet), pointer:: head => null()
        integer(kind=4), allocatable, dimension(:):: current_j

        this%nnz = src%nnz
        this%n_rows = src%n_rows
        this%n_cols = src%n_cols

        allocate(this%rows(this%n_rows+1))
        allocate(this%cols(this%nnz))
        allocate(this%vals(this%nnz))
        allocate(current_j(this%n_rows))

        ! note that Fortran is 1-based, not 0-based, so it's a little tricky
        this%rows(1) = 1

        ! fisrt loop to get nnz in each row
        head => src%bg
        do i_nnz = 1, this%nnz
            i = head%idx(1)
            this%rows(i+1) = this%rows(i+1) + 1
            head => head%next
        end do
        nullify(head)

        ! accumulate
        do i = 2, this%n_rows + 1
            this%rows(i) = this%rows(i) + this%rows(i-1)
        end do

        ! check
        if (this%rows(this%n_rows+1) /= (this%nnz + 1)) then
            print *, "Error in init_CSR."
            stop
        end if

        ! copy the current state of column index counter
        current_j = this%rows(1:this%n_rows)

        ! copy column indices and values to CSR
        head => src%bg
        do i_nnz = 1, this%nnz
            i = head%idx(1)
            this%cols(current_j(i)) = head%idx(2)
            current_j(i) = current_j(i) + 1
            head => head%next
        end do
        nullify(head)

        ! check
        do i = 1, this%n_rows
            if (current_j(i) /= this%rows(i+1)) then
                print *, "Error in init_CSR: current_j(", i, ")=", &
                    current_j(i), "!= this%rows(", i+1, ")=", this%rows(i+1)
                stop
            end if
        end do

        ! free memory
        deallocate(current_j)

        ! sort columns
        do i = 1, this%n_rows
            call qsort(this%cols(this%rows(i):this%rows(i+1)-1), &
                int(this%rows(i+1)-this%rows(i), c_size_t), &
                c_sizeof(this%cols(i)), compare_int4)
        end do

        ! set values according to sorted columns
        head => src%bg
        do i_nnz = 1, this%nnz
            i = head%idx(1)
            where(this%cols(this%rows(i):this%rows(i+1)-1) == head%idx(2))
                this%vals(this%rows(i):this%rows(i+1)-1) = head%val
            end where
            head => head%next
        end do
        nullify(head)
    end subroutine init_CSR

    ! get_CSR
    function get_CSR(this, i, j) result(val)
        use, intrinsic:: iso_c_binding, only: c_ptr, c_size_t, c_associated
        use, intrinsic:: iso_c_binding, only: c_f_pointer, c_sizeof
        class(CSR), intent(in):: this
        integer(kind=4), intent(in):: i, j
        real(kind=8):: val

        ! local variables
        type(c_ptr):: c_loc_ptr
        integer(kind=4), pointer:: f_loc_ptr
        integer(kind=4):: loc_cols

        c_loc_ptr = bsearch(j, this%cols(this%rows(i):this%rows(i+1)-1), &
            int(this%rows(i+1)-this%rows(i), c_size_t), &
            c_sizeof(this%cols(this%rows(i))), compare_int4)

        if (.not. c_associated(c_loc_ptr)) then
            val = 0D0
            return
        end if

        call c_f_pointer(c_loc_ptr, f_loc_ptr)
        loc_cols = (loc(f_loc_ptr) - loc(this%cols(1))) / sizeof(this%cols(1)) + 1
        val = this%vals(loc_cols)
    end function get_CSR

    ! set_CSR
    subroutine set_CSR(this, i, j, val)
        use, intrinsic:: iso_c_binding, only: c_ptr, c_size_t, c_associated
        use, intrinsic:: iso_c_binding, only: c_f_pointer, c_sizeof
        class(CSR), intent(inout):: this
        integer(kind=4), intent(in):: i, j
        real(kind=8), intent(in):: val

        ! local variables
        type(c_ptr):: c_loc_ptr
        integer(kind=4), pointer:: f_loc_ptr
        integer(kind=4):: loc_cols

        c_loc_ptr = bsearch(j, this%cols(this%rows(i):this%rows(i+1)-1), &
            int(this%rows(i+1)-this%rows(i), c_size_t), &
            c_sizeof(this%cols(this%rows(i))), compare_int4)

        if (.not. c_associated(c_loc_ptr)) then
            print *, "Error: trying to set value to zero location in CSR."
            stop
        end if

        call c_f_pointer(c_loc_ptr, f_loc_ptr)
        loc_cols = (loc(f_loc_ptr) - loc(this%cols(1))) / sizeof(this%cols(1)) + 1
        this%vals(loc_cols) = val
    end subroutine set_CSR

    ! add_CSR
    subroutine add_CSR(this, i, j, val)
        use, intrinsic:: iso_c_binding, only: c_ptr, c_size_t, c_associated
        use, intrinsic:: iso_c_binding, only: c_f_pointer, c_sizeof
        class(CSR), intent(inout):: this
        integer(kind=4), intent(in):: i, j
        real(kind=8), intent(in):: val

        ! local variables
        type(c_ptr):: c_loc_ptr
        integer(kind=4), pointer:: f_loc_ptr
        integer(kind=4):: loc_cols

        c_loc_ptr = bsearch(j, this%cols(this%rows(i):this%rows(i+1)-1), &
            int(this%rows(i+1)-this%rows(i), c_size_t), &
            c_sizeof(this%cols(this%rows(i))), compare_int4)

        if (.not. c_associated(c_loc_ptr)) then
            print *, "Error: trying to set value to zero location in CSR."
            stop
        end if

        call c_f_pointer(c_loc_ptr, f_loc_ptr)
        loc_cols = (loc(f_loc_ptr) - loc(this%cols(1))) / sizeof(this%cols(1)) + 1
        this%vals(loc_cols) = this%vals(loc_cols) + val
    end subroutine add_CSR

    ! count_block_CSR
    function count_block_CSR(this, rowl, rowh, coll, colh) result(ans)
        class(CSR), intent(in):: this
        integer(kind=4), intent(in):: rowl, rowh, coll, colh
        integer(kind=4):: ans

        integer(kind=4):: i
        integer(kind=4):: low_bound, high_bound

        ans = 0
        do i = rowl, rowh
            low_bound = this%rows(i)
            do while (low_bound<this%rows(i+1))
                if (this%cols(low_bound) >= coll) exit
                low_bound = low_bound + 1
            end do

            ! all col indices in this row are smaller than coll
            if (low_bound == this%rows(i+1)) exit ! go to next row
            
            high_bound = low_bound
            do while (high_bound<this%rows(i+1))
                if (this%cols(high_bound) <= colh) then
                    ans = ans + 1
                    high_bound = high_bound + 1
                else
                    exit
                end if
            end do
        end do
    end function count_block_CSR

    ! count_block_CSR
    subroutine add_block_CSR(this, rowl, rowh, coll, colh, val)
        class(CSR), intent(inout):: this
        integer(kind=4), intent(in):: rowl, rowh, coll, colh
        real(kind=8), intent(in):: val

        integer(kind=4):: i
        integer(kind=4):: low_bound, high_bound

        do i = rowl, rowh
            low_bound = this%rows(i)
            do while (low_bound<this%rows(i+1))
                if (this%cols(low_bound) >= coll) exit
                low_bound = low_bound + 1
            end do

            ! all col indices in this row are smaller than coll
            if (low_bound == this%rows(i+1)) exit ! go to next row
            
            high_bound = low_bound
            do while (high_bound<this%rows(i+1))
                if (this%cols(high_bound) <= colh) then
                    this%vals(high_bound) = this%vals(high_bound) + val
                    high_bound = high_bound + 1
                else
                    exit
                end if
            end do
        end do
    end subroutine add_block_CSR

    ! write_CSR
    subroutine write_CSR(this, iounit, iotype, v_list, stat, msg)
        ! variable declaration
        class(CSR), intent(in):: this
        integer(kind=4), intent(in):: iounit
        character(*), intent(in)::iotype
        integer(kind=4), intent(in):: v_list(:)
        integer(kind=4), intent(out):: stat
        character(*), intent(inout):: msg
        integer(kind=4):: i, j
        character:: n, t


        n = new_line(t) ! n is the "new line" character
        t = achar(9) ! t is the character for a tab

        write(iounit, *, iostat=stat, iomsg=msg) "NNZ = ", this%nnz
        write(iounit, *, iostat=stat, iomsg=msg) n

        do i = 1, this%n_rows
            do j = this%rows(i), this%rows(i+1)-1
                write(iounit, *, iostat=stat, iomsg=msg) &
                    "(", i, ",", this%cols(j), ") =>", this%vals(j)
                write(iounit, *, iostat=stat, iomsg=msg) n
            end do
        end do
    end subroutine write_CSR

    ! destroy_CSR
    subroutine destroy_CSR(this)
        class(CSR), intent(inout):: this
        type(IndexSet), pointer:: head

        this%nnz = 0
        this%n_rows = 0
        this%n_cols = 0

        if (allocated(this%rows)) deallocate(this%rows)
        if (allocated(this%cols)) deallocate(this%cols)
        if (allocated(this%vals)) deallocate(this%vals)
    end subroutine destroy_CSR

    ! destructor_CSR
    subroutine destructor_CSR(this)
        type(CSR), intent(inout):: this
        call this%destroy()
    end subroutine destructor_CSR

    ! compare_int4
    function compare_int4(a, b)
        use, intrinsic:: iso_c_binding, only: c_int, c_ptr, c_f_pointer
        type(c_ptr), intent(in), value:: a, b
        integer(kind=c_int):: compare_int4
        integer(kind=c_int), pointer:: ap, bp

        call c_f_pointer(a, ap)
        call c_f_pointer(b, bp)

        compare_int4 = int(ap-bp, c_int)
    end function compare_int4

    ! compress
    module subroutine compress(A, sp)
        real(kind=8), dimension(:, :), intent(in):: A
        type(CSR), intent(inout):: sp

        ! local variables
        integer(kind=4):: i, j

        ! clear sp
        call sp%destroy()

        sp%nnz = count(dabs(A) >= 1e-9)
        sp%n_rows = size(A, 1)
        sp%n_cols = size(A, 2)

        allocate(sp%rows(sp%n_rows+1))
        allocate(sp%cols(sp%nnz))
        allocate(sp%vals(sp%nnz))

        sp%rows(1) = 1

        ! TODO: bad implementation. Expensive. Improve in the future.
        do i = 1, sp%n_rows
            sp%rows(i+1) = count(dabs(A(i, :))>=1e-9) + sp%rows(i)

            sp%cols(sp%rows(i):sp%rows(i+1)-1) = &
                pack((/ (j, j=1, sp%n_cols) /), dabs(A(i, :))>=1e-9)

            sp%vals(sp%rows(i):sp%rows(i+1)-1) = &
                pack(A(i, :), dabs(A(i, :))>=1e-9)
        end do

    end subroutine compress

end module SPM_module
