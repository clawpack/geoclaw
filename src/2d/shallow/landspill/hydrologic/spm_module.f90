!
! SPM_module.f90
! Contributors: Pi-Yueh Chuang <pychuang@gwu.edu>
!
! Distributed under terms of the BSD 3-Clause license.
!

!> @brief Naive sparse-matrix implementations.
module SPM_module
    implicit none
    private
    public:: IndexSet, COO, MultiLayerCSR, compress

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

    !> @brief A naive implementation of the sparse pattern of CSR format.
    type:: SparsePattern
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

        contains
        !> @brief Init the sparse pattern from a COO matrix.
        procedure:: init => init_SparsePattern
        !> @brief Given a (i, j), return the true 1d index in CSR format.
        procedure:: get_index => get_index_SparsePattern
        !> @brief Given multiple (i, j)s, return the true 1d indices in CSR format.
        procedure:: get_indices => get_indices_SparsePattern
        !> @brief Destroy this SparsePattern.
        procedure:: destroy => destroy_SparsePattern
        !> @brief Write.
        procedure, private:: write_SparsePattern
        !> @brief Overloading intrinsic write.
        generic:: write(formatted) => write_SparsePattern
        !> @brief Destructor.
        final:: destructor_SparsePattern
    end type SparsePattern

    !> @brief A naive implementation of multi-layer CSR sparse matrix.
    type:: MultiLayerCSR
        private
        !> @brief Underlying sparse pattern.
        type(SparsePattern):: sp
        !> @brief Number of layers in the matrix.
        integer(kind=4):: n_layers = 0
        !> @brief Values of non-zeros.
        real(kind=8), allocatable, dimension(:, :):: vals

        contains
        !> @brief Init from a MultiLayerCSR matrix.
        procedure:: init => init_MultiLayerCSR
        !> @brief Given a (i, j), return the true 1d index in CSR format.
        procedure:: get_index => get_index_MultiLayerCSR
        !> @brief Given multiple (i, j)s, return the true 1d indices in CSR format.
        procedure:: get_indices => get_indices_MultiLayerCSR
        !> @brief Access the value of a matrix element.
        procedure:: get_value => get_value_MultiLayerCSR
        !> @brief Access the values of multiple matrix elements.
        procedure:: get_values => get_values_MultiLayerCSR
        !> @brief Set a matrix element.
        procedure:: set => set_MultiLayerCSR
        !> @brief Set multiple matrix elements to the same value.
        procedure:: set_multiples => set_multiples_MultiLayerCSR
        !> @brief Set a single value to all non-zero elements.
        procedure:: set_all => set_all_MultiLayerCSR
        !> @brief Add a value to a matrix element.
        procedure:: add => add_MultiLayerCSR
        !> @brief Add a single value to multiple matrix elements.
        procedure:: add_multiples => add_multiples_MultiLayerCSR
        !> @brief Add a single value to all non-zero elements.
        procedure:: add_all => add_all_MultiLayerCSR
        !> @brief In-place multiplying a value to a matrix element.
        procedure:: mult => mult_MultiLayerCSR
        !> @brief In-place multiplying a value to multiple matrix elements.
        procedure:: mult_multiples => mult_multiples_MultiLayerCSR
        !> @brief In-place multiplying a value to all non-zero elements.
        procedure:: mult_all => mult_all_MultiLayerCSR
        !> @brief Summation.
        procedure:: sum => sum_MultiLayerCSR
        !> @brief Destroy this MultiLayerCSR.
        procedure:: destroy => destroy_MultiLayerCSR
        !> @brief Write.
        procedure, private:: write_MultiLayerCSR
        !> @brief Overloading intrinsic write.
        generic:: write(formatted) => write_MultiLayerCSR
        !> @brief Destructor.
        final:: destructor_MultiLayerCSR
    end type MultiLayerCSR

    !> @brief An interface to standard C qsort function
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

    !> @brief An interface to standard C bsearch function
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

    !> @brief An explicit interface so that Intel compiler can work.
    interface compress
        module procedure:: compress
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

    ! constructor_SparsePattern
    function constructor_SparsePattern(src) result(sp)
        type(COO), intent(in):: src
        type(SparsePattern):: sp

        call init_SparsePattern(sp, src)
    end function constructor_SparsePattern

    ! init_SparsePattern
    subroutine init_SparsePattern(this, src)
        use, intrinsic:: iso_c_binding, only: c_size_t, c_sizeof
        class(SparsePattern), intent(inout):: this
        type(COO), intent(in):: src
        integer(kind=4):: i, j, i_nnz
        type(IndexSet), pointer:: head => null()
        integer(kind=4), allocatable, dimension(:):: current_j

        this%nnz = src%nnz
        this%n_rows = src%n_rows
        this%n_cols = src%n_cols

        allocate(this%rows(this%n_rows+1))
        allocate(this%cols(this%nnz))
        allocate(current_j(this%n_rows))

        ! note that Fortran is 1-based, not 0-based, so it's a little tricky
        this%rows = 0
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

        ! copy column indices to CSR
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
    end subroutine init_SparsePattern

    ! get_index_SparsePattern
    function get_index_SparsePattern(this, i, j) result(idx)
        use, intrinsic:: iso_c_binding, only: c_ptr, c_size_t, c_associated
        use, intrinsic:: iso_c_binding, only: c_f_pointer, c_sizeof
        class(SparsePattern), intent(in):: this
        integer(kind=4), intent(in):: i, j
        integer(kind=4):: idx

        ! local variables
        type(c_ptr):: c_loc_ptr
        integer(kind=4), pointer:: f_loc_ptr

        c_loc_ptr = bsearch(j, this%cols(this%rows(i):this%rows(i+1)-1), &
            int(this%rows(i+1)-this%rows(i), c_size_t), &
            c_sizeof(this%cols(this%rows(i))), compare_int4)

        if (.not. c_associated(c_loc_ptr)) then
            idx = 0D0
            return
        end if

        call c_f_pointer(c_loc_ptr, f_loc_ptr)
        idx = (loc(f_loc_ptr) - loc(this%cols(1))) / sizeof(this%cols(1)) + 1
    end function get_index_SparsePattern

    ! get_indices_SparsePattern
    subroutine get_indices_SparsePattern(this, rowl, rowh, coll, colh, nnzs, idxs)
        class(SparsePattern), intent(in):: this
        integer(kind=4), intent(in):: rowl, rowh, coll, colh
        integer(kind=4), intent(out):: nnzs
        integer(kind=4), dimension(:), intent(out):: idxs

        integer(kind=4):: i
        integer(kind=4):: low_bound, high_bound

        nnzs = 0
        idxs = 0
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
                    nnzs = nnzs + 1
                    idxs(nnzs) = high_bound
                    high_bound = high_bound + 1
                else
                    exit
                end if
            end do
        end do
    end subroutine get_indices_SparsePattern

    ! write_SparsePattern
    subroutine write_SparsePattern(this, iounit, iotype, v_list, stat, msg)
        ! variable declaration
        class(SparsePattern), intent(in):: this
        integer(kind=4), intent(in):: iounit
        character(*), intent(in)::iotype
        integer(kind=4), intent(in):: v_list(:)
        integer(kind=4), intent(out):: stat
        character(*), intent(inout):: msg
        integer(kind=4):: i, j, k
        character:: n, t


        n = new_line(t) ! n is the "new line" character
        t = achar(9) ! t is the character for a tab

        write(iounit, *, iostat=stat, iomsg=msg) "NNZ = ", this%nnz
        write(iounit, *, iostat=stat, iomsg=msg) n

        do i = 1, this%n_rows
            do j = this%rows(i), this%rows(i+1)-1
                write(iounit, *, iostat=stat, iomsg=msg) &
                    "(", i, ",", this%cols(j), ")"
                write(iounit, *, iostat=stat, iomsg=msg) n
            end do
        end do
    end subroutine write_SparsePattern

    ! destroy_SparsePattern
    subroutine destroy_SparsePattern(this)
        class(SparsePattern), intent(inout):: this

        this%nnz = 0
        this%n_rows = 0
        this%n_cols = 0

        if (allocated(this%rows)) deallocate(this%rows)
        if (allocated(this%cols)) deallocate(this%cols)
    end subroutine destroy_SparsePattern

    ! destructor_SparsePattern
    subroutine destructor_SparsePattern(this)
        type(SparsePattern), intent(inout):: this
        call this%destroy()
    end subroutine destructor_SparsePattern

    ! constructor_MultiLayerCSR
    function constructor_MultiLayerCSR(n_layers, src) result(mtx)
        type(COO), intent(in):: src
        integer(kind=4), intent(in):: n_layers
        type(MultiLayerCSR):: mtx

        call init_MultiLayerCSR(mtx, n_layers, src)
    end function constructor_MultiLayerCSR

    ! init_MultiLayerCSR
    subroutine init_MultiLayerCSR(this, n_layers, src)
        use, intrinsic:: iso_c_binding, only: c_size_t, c_sizeof
        class(MultiLayerCSR), intent(inout):: this
        integer(kind=4), intent(in):: n_layers
        type(COO), intent(in):: src
        integer(kind=4):: i, i_nnz
        type(IndexSet), pointer:: head => null()

        call this%sp%init(src)
        this%n_layers = n_layers

        allocate(this%vals(this%sp%nnz, this%n_layers))
        this%vals = 0D0

        ! set values according to sorted columns
        head => src%bg
        do i_nnz = 1, this%sp%nnz
            i = head%idx(1)
            where(this%sp%cols(this%sp%rows(i):this%sp%rows(i+1)-1) == head%idx(2))
                this%vals(this%sp%rows(i):this%sp%rows(i+1)-1, 1) = head%val
            end where
            head => head%next
        end do
        nullify(head)
    end subroutine init_MultiLayerCSR

    ! get_index_MultiLayerCSR
    function get_index_MultiLayerCSR(this, i, j) result(idx)
        class(MultiLayerCSR), intent(in):: this
        integer(kind=4), intent(in):: i, j
        integer(kind=4):: idx

        idx = this%sp%get_index(i, j)
    end function get_index_MultiLayerCSR

    ! get_indices_MultiLayerCSR
    subroutine get_indices_MultiLayerCSR(this, rowl, rowh, coll, colh, nnzs, idxs)
        class(MultiLayerCSR), intent(in):: this
        integer(kind=4), intent(in):: rowl, rowh, coll, colh
        integer(kind=4), intent(out):: nnzs
        integer(kind=4), dimension(:), intent(out):: idxs

        call this%sp%get_indices(rowl, rowh, coll, colh, nnzs, idxs)
    end subroutine get_indices_MultiLayerCSR

    ! get_value_MultiLayerCSR
    function get_value_MultiLayerCSR(this, idx, layer) result(val)
        class(MultiLayerCSR), intent(in):: this
        integer(kind=4), intent(in):: idx, layer
        real(kind=8):: val

        val = this%vals(idx, layer)
    end function get_value_MultiLayerCSR

    ! get_values_MultiLayerCSR
    subroutine get_values_MultiLayerCSR(this, nnzs, idxs, layer, vals)
        class(MultiLayerCSR), intent(in):: this
        integer(kind=4), intent(in):: nnzs, layer
        integer(kind=4), dimension(nnzs), intent(in):: idxs
        real(kind=8), dimension(nnzs), intent(out):: vals

        vals = this%vals(idxs, layer)
    end subroutine get_values_MultiLayerCSR

    ! set_MultiLayerCSR
    subroutine set_MultiLayerCSR(this, idx, layer, val)
        class(MultiLayerCSR), intent(inout):: this
        integer(kind=4), intent(in):: idx, layer
        real(kind=8), intent(in):: val

        this%vals(idx, layer) = val
    end subroutine set_MultiLayerCSR

    ! set_multiples_MultiLayerCSR
    subroutine set_multiples_MultiLayerCSR(this, nnzs, idxs, layer, val)
        class(MultiLayerCSR), intent(inout):: this
        integer(kind=4), intent(in):: nnzs, layer
        integer(kind=4), dimension(nnzs), intent(in):: idxs
        real(kind=8), intent(in):: val

        this%vals(idxs, layer) = val
    end subroutine set_multiples_MultiLayerCSR

    ! set_all_MultiLayerCSR
    subroutine set_all_MultiLayerCSR(this, layer, val)
        class(MultiLayerCSR), intent(inout):: this
        integer(kind=4), intent(in):: layer
        real(kind=8), intent(in):: val

        this%vals(:, layer) = val
    end subroutine set_all_MultiLayerCSR

    ! add_MultiLayerCSR
    subroutine add_MultiLayerCSR(this, idx, layer, val)
        class(MultiLayerCSR), intent(inout):: this
        integer(kind=4), intent(in):: idx, layer
        real(kind=8), intent(in):: val

        this%vals(idx, layer) = this%vals(idx, layer) + val
    end subroutine add_MultiLayerCSR

    ! add_multiples_MultiLayerCSR
    subroutine add_multiples_MultiLayerCSR(this, nnzs, idxs, layer, val)
        class(MultiLayerCSR), intent(inout):: this
        integer(kind=4), intent(in):: nnzs, layer
        integer(kind=4), dimension(nnzs), intent(in):: idxs
        real(kind=8), intent(in):: val

        this%vals(idxs, layer) = this%vals(idxs, layer) + val
    end subroutine add_multiples_MultiLayerCSR

    ! add_all_MultiLayerCSR
    subroutine add_all_MultiLayerCSR(this, layer, val)
        class(MultiLayerCSR), intent(inout):: this
        integer(kind=4), intent(in):: layer
        real(kind=8), intent(in):: val

        this%vals(:, layer) = this%vals(:, layer) + val
    end subroutine add_all_MultiLayerCSR

    ! mult_MultiLayerCSR
    subroutine mult_MultiLayerCSR(this, idx, layer, val)
        class(MultiLayerCSR), intent(inout):: this
        integer(kind=4), intent(in):: idx, layer
        real(kind=8), intent(in):: val

        this%vals(idx, layer) = this%vals(idx, layer) * val
    end subroutine mult_MultiLayerCSR

    ! mult_multiples_MultiLayerCSR
    subroutine mult_multiples_MultiLayerCSR(this, nnzs, idxs, layer, val)
        class(MultiLayerCSR), intent(inout):: this
        integer(kind=4), intent(in):: nnzs, layer
        integer(kind=4), dimension(nnzs), intent(in):: idxs
        real(kind=8), intent(in):: val

        this%vals(idxs, layer) = this%vals(idxs, layer) * val
    end subroutine mult_multiples_MultiLayerCSR

    ! mult_all_MultiLayerCSR
    subroutine mult_all_MultiLayerCSR(this, layer, val)
        class(MultiLayerCSR), intent(inout):: this
        integer(kind=4), intent(in):: layer
        real(kind=8), intent(in):: val

        this%vals(:, layer) = this%vals(:, layer) * val
    end subroutine mult_all_MultiLayerCSR

    ! sum_MultiLayerCSR
    function sum_MultiLayerCSR(this, layer) result(ans)
        class(MultiLayerCSR), intent(in):: this
        integer(kind=4), intent(in):: layer
        real(kind=8):: ans

        ans = sum(this%vals(:, layer))
    end function sum_MultiLayerCSR

    ! write_MultiLayerCSR
    subroutine write_MultiLayerCSR(this, iounit, iotype, v_list, stat, msg)
        ! variable declaration
        class(MultiLayerCSR), intent(in):: this
        integer(kind=4), intent(in):: iounit
        character(*), intent(in)::iotype
        integer(kind=4), intent(in):: v_list(:)
        integer(kind=4), intent(out):: stat
        character(*), intent(inout):: msg
        integer(kind=4):: i, j, k
        character:: n, t


        n = new_line(t) ! n is the "new line" character
        t = achar(9) ! t is the character for a tab

        write(iounit, *, iostat=stat, iomsg=msg) "NNZ = ", this%sp%nnz
        write(iounit, *, iostat=stat, iomsg=msg) "N Layers = ", this%n_layers
        write(iounit, *, iostat=stat, iomsg=msg) n

        do k = 1, this%n_layers
            write(iounit, *, iostat=stat, iomsg=msg) "Layer", k
            do i = 1, this%sp%n_rows
                do j = this%sp%rows(i), this%sp%rows(i+1)-1
                    write(iounit, *, iostat=stat, iomsg=msg) &
                        "(", i, ",", this%sp%cols(j), ") =>", this%vals(j, k)
                    write(iounit, *, iostat=stat, iomsg=msg) n
                end do
            end do
        end do
    end subroutine write_MultiLayerCSR

    ! destroy_MultiLayerCSR
    subroutine destroy_MultiLayerCSR(this)
        class(MultiLayerCSR), intent(inout):: this

        this%n_layers = 0
        if (allocated(this%vals)) deallocate(this%vals)
        call this%sp%destroy()
    end subroutine destroy_MultiLayerCSR

    ! destructor_MultiLayerCSR
    subroutine destructor_MultiLayerCSR(this)
        type(MultiLayerCSR), intent(inout):: this
        call this%destroy()
    end subroutine destructor_MultiLayerCSR

    ! compare_int4
    function compare_int4(a, b) bind(C)
        use, intrinsic:: iso_c_binding, only: c_int, c_ptr, c_f_pointer
        type(c_ptr), intent(in), value:: a, b
        integer(kind=c_int):: compare_int4
        integer(kind=c_int), pointer:: ap, bp

        call c_f_pointer(a, ap)
        call c_f_pointer(b, bp)

        compare_int4 = int(ap-bp, c_int)
    end function compare_int4

    ! compress
    subroutine compress(A, mtx)
        real(kind=8), dimension(:, :), intent(in):: A
        type(MultiLayerCSR), intent(inout):: mtx

        ! local variables
        integer(kind=4):: i, j

        ! clear sp
        call mtx%destroy()

        mtx%sp%nnz = count(dabs(A) >= 1e-9)
        mtx%sp%n_rows = size(A, 1)
        mtx%sp%n_cols = size(A, 2)
        mtx%n_layers = 1

        allocate(mtx%sp%rows(mtx%sp%n_rows+1))
        allocate(mtx%sp%cols(mtx%sp%nnz))
        allocate(mtx%vals(mtx%sp%nnz, mtx%n_layers))

        mtx%sp%rows(1) = 1

        ! TODO: bad implementation. Expensive. Improve in the future.
        do i = 1, mtx%sp%n_rows
            mtx%sp%rows(i+1) = count(dabs(A(i, :))>=1e-9) + mtx%sp%rows(i)

            mtx%sp%cols(mtx%sp%rows(i):mtx%sp%rows(i+1)-1) = &
                pack((/ (j, j=1, mtx%sp%n_cols) /), dabs(A(i, :))>=1e-9)

            mtx%vals(mtx%sp%rows(i):mtx%sp%rows(i+1)-1, 1) = &
                pack(A(i, :), dabs(A(i, :))>=1e-9)
        end do

    end subroutine compress

end module SPM_module
