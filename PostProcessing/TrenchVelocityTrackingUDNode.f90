program trench_velocity_tracking_udnode
  implicit none

  integer, parameter :: max_elements = 136000
  integer, parameter :: max_nodes    = 400000

  integer :: v, p, f, i, j, stat, pos
  integer :: elem_id, mat, node_id
  real(8) :: x, y, time, max_x_UCC
  integer :: id_te

  integer :: conn(max_elements,4)
  real(8) :: vx, vy

  character(len=100)  :: filename, nodi_file
  character(len=100)  :: output_filename_OC
  character(len=1000) :: base_dir, grid_dir, tracking_dir
  character(len=200)  :: line
  logical :: exists

  character(len=20) :: viscosity(3), plasticity(3), forcing(25)
  character(len=8)  :: visco_short(3), plast_short(3), force_short(25)

  real(8) :: a,b,c,d,e,z,g,h,o,u,k,l,m,n

  viscosity   = (/'B_ev5-10 ', 'C_ev02-05', 'D_ev05-15'/)
  visco_short = (/'B', 'C', 'D'/)
  plasticity   = (/'Peierls_PW2_VSall', 'Peierls_PW3_VSall', 'Peierls_PW4_VSall'/)
  plast_short  = (/'PW2', 'PW3', 'PW4'/)
  forcing = (/'Spontaneous     ','Forced/v1       ','Forced/v05      ','Forced/v025     ', &
              'Forced/v01      ','Forced/v005     ','Forced/v001     ','PostSoft10/v1   ', &
              'PostSoft10/v05  ','PostSoft10/v025 ','PostSoft10/v01  ','PostSoft10/v005 ', &
              'PostSoft10/v001 ','PostSoft20/v1   ','PostSoft20/v05  ','PostSoft20/v025 ', &
              'PostSoft20/v01  ','PostSoft20/v005 ','PostSoft20/v001 ','PostSoft30/v1   ', &
              'PostSoft30/v05  ','PostSoft30/v025 ','PostSoft30/v01  ','PostSoft30/v005 ', &
              'PostSoft30/v001 '/)
  force_short = (/'spont   ','f_v1    ','f_v05   ','f_v025  ', &
                  'f_v01   ','f_v005  ','f_v001  ','p1_v1   ','p1_v05  ', &
                  'p1_v025 ','p1_v01  ','p1_v005 ','p1_v001 ','p2_v1   ','p2_v05  ', &
                  'p2_v025 ','p2_v01  ','p2_v005 ','p2_v001 ','p3_v1   ','p3_v05  ', &
                  'p3_v025 ','p3_v01  ','p3_v005 ','p3_v001 '/)

  base_dir     = "/media/valeria/LaCie/BACKUP_PHD_FEDELI/DOTTORATO/SI_passiveMargin/Output_models/VF_PassiveMargin/"
  tracking_dir = "Tracking/"

  open(unit=99, file=trim(base_dir)//"Griglia_final", status='old', iostat=stat)
  if (stat /= 0) then
    print *, "Errore: impossibile aprire Griglia_final"
    stop
  end if
  do i = 1, 136977
    read(99,*)
  end do
  do i = 1, max_elements
    read(99,*) elem_id, conn(i,1), conn(i,2), conn(i,3), conn(i,4)
  end do
  close(99)

  do v = 1, 3
    do p = 1, 3
      do f = 1, 25
        grid_dir = trim(base_dir)//trim(viscosity(v))//"/20Ma_Eq3/"//trim(plasticity(p))//"/"// &
                   trim(forcing(f))//"/txt_files/Grid/"

        inquire(file=trim(grid_dir)//"Elementi.00000.txt", exist=exists)
        if (.not. exists) then
          print *, "Cartella non trovata: ", visco_short(v), plast_short(p), force_short(f)
          cycle
        end if

        output_filename_OC = trim(tracking_dir)//"Trench_"//trim(visco_short(v))//"_"//trim(plast_short(p))&
                              //"_"//trim(force_short(f))//".txt"
        inquire(file=output_filename_OC, exist=exists)
        if (exists) then
          print*, "File esiste già: ",  visco_short(v), plast_short(p), force_short(f)
          cycle
        endif

        print*, "> ",  visco_short(v), plast_short(p), force_short(f)
        open(unit=13,file=output_filename_OC ,status='unknown')
        write(13,*) 'Time [Myr], X [m], Y [m], Vx [cm/yr], Vy [cm/yr], V [cm/yr]'

        do i = 0, 60            
          write(filename,'(A,I5.5,A)') 'Elementi.', i, '.txt'
          open(unit=20,file=trim(grid_dir)//trim(filename),status='old',iostat=stat)
          if (stat/=0) exit

          read(20,'(A)') line
          pos = scan(line,'0123456789')
          if (pos>0) read(line(pos:),*) time
          do j=1,22; read(20,*); end do

          max_x_UCC = 0.0D0

          do
            read(20,*,iostat=stat) elem_id, x, y, a,b,c,d,e,z,g,h,o,u,k,l,m,n,mat
            if (stat/=0) exit
            if ((mat == 7) .and. (x > max_x_UCC)) then
                max_x_UCC = x
                id_te = elem_id
            end if
          end do
          close(20)

          !print*, id_te, max_x_UCC, conn(id_te,3)
          write(nodi_file,'(A,I5.5,A)') 'Nodi.', i, '.txt'
          open(unit=30,file=trim(grid_dir)//trim(nodi_file),status='old',iostat=stat)
          if (stat/=0) cycle
          do j=1,10; read(30,*); end do
          

          do            
            read(30,*,iostat=stat) node_id, x, y, vx, vy
            if (stat/=0) exit
            !print*, time, node_id, x, y, vx!, vy
            if (node_id == conn(id_te,3)) then
              !PRINT*, node_id, time, x, y, vx, vy, sqrt(vx**2 + vy**2)
              write(13, '(F10.3,2F12.3,3F12.3)') time, x, y, vx, vy, sqrt(vx**2 + vy**2)
            end if
            
          end do
          close(30)

        end do

        close(13)

      end do
    end do
  end do

end program 
