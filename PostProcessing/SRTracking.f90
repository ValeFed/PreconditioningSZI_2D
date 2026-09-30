program sr_tracking
  implicit none
  integer, parameter :: max_nodes = 136971, max_elements = 136000
  integer :: i, j, node, stat, pos, v, p, f, inode, elem_id
  integer :: conn(max_elements,4)
  real, dimension(max_nodes) :: xn, yn
  real(8) :: x, y, time, Nel, Nemin, sr, xT, tT
  real(8) :: srmin,a4,a5,a6,a7,a8, SRmean
  real(8), parameter :: xmin = 350000.D0, xmax = 10000.D0, ymin = 510000.0D0, ymax = 710000.0D0
  character(len=100) :: filename, output_filename, trench_filename
  character(len=1000) :: base_dir, grid_dir, tracking_dir
  character(len=200) :: line
  character(len=20) :: viscosity(3), plasticity(3), forcing(25)
  character(len=8)  :: visco_short(3), plast_short(3), force_short(25)
  logical :: file_exists, dir_exists

  ! Definizione delle possibili combinazioni
  viscosity = (/'B_ev5-10 ', 'C_ev02-05', 'D_ev05-15'/)
  visco_short = (/'B', 'C', 'D'/)

  plasticity = (/'Peierls_PW2_VSall', 'Peierls_PW3_VSall', 'Peierls_PW4_VSall'/)
  plast_short = (/'PW2', 'PW3', 'PW4'/)

  srmin=2.e-14
  
  forcing = (/'Spontaneous     ', 'Forced/v1       ', 'Forced/v05      ', 'Forced/v025     ', &
  'Forced/v01      ', 'Forced/v005     ', 'Forced/v001     ', 'PostSoft10/v1   ',  &
  'PostSoft10/v05  ', 'PostSoft10/v025 ', 'PostSoft10/v01  ', 'PostSoft10/v005 ', &
  'PostSoft10/v001 ', 'PostSoft20/v1   ', 'PostSoft20/v05  ', 'PostSoft20/v025 ', &
  'PostSoft20/v01  ', 'PostSoft20/v005 ', 'PostSoft20/v001 ', 'PostSoft30/v1   ',  &                
  'PostSoft30/v05  ', 'PostSoft30/v025 ', 'PostSoft30/v01  ', 'PostSoft30/v005 ', &
  'PostSoft30/v001 '/)
force_short = (/'spont   ', 'f_v1    ', 'f_v05   ', 'f_v025  ', &
      'f_v01   ', 'f_v005  ', 'f_v001  ', 'p1_v1   ', 'p1_v05  ',&
      'p1_v025 ', 'p1_v01  ', 'p1_v005 ', 'p1_v001 ', 'p2_v1   ', 'p2_v05  ',&
      'p2_v025 ', 'p2_v01  ', 'p2_v005 ', 'p2_v001 ', 'p3_v1   ', 'p3_v05  ',&
      'p3_v025 ', 'p3_v01  ', 'p3_v005 ', 'p3_v001 '/)

  ! Definizione delle cartelle
  base_dir = "/media/valeria/LaCie/BACKUP_PHD_FEDELI/DOTTORATO/SI_passiveMargin/Output_models/VF_PassiveMargin/"
  tracking_dir = "Tracking/"

  open(unit=99, file=trim(base_dir)//"Griglia_final", status='old', iostat=stat)
  if (stat /= 0) then
    print *, "Errore: impossibile aprire Griglia_final"
    stop
  end if

  do i=1,5; read(99,*); end do  ! Salta le prime 5 righe
  do i = 1, max_nodes
    read(99,*) inode, xn(inode), yn(inode)
  end do
  read(99,*)
  do i = 1, max_elements
    read(99,*) elem_id, conn(i,1), conn(i,2), conn(i,3), conn(i,4)
  end do
  close(99)

  ! Loop sulle combinazioni di modelli
  do v = 1, 3
    do p = 1, 3
      do f = 1, 25
        ! Creazione del percorso della cartella Grid
        grid_dir = trim(base_dir) // trim(viscosity(v)) // "/20Ma_Eq3/" // trim(plasticity(p)) // "/" // &
                   trim(forcing(f)) // "/txt_files/Grid/"
        trench_filename = trim(tracking_dir)//"Trench_"//trim(visco_short(v))//"_"//trim(plast_short(p))&
                              //"_"//trim(force_short(f))//".txt"
        ! Verifica se la cartella esiste
        inquire(file=grid_dir, exist=dir_exists)
        if (.not. dir_exists) then
          print *, "La cartella non esiste: ", visco_short(v), plast_short(p), force_short(f)
          cycle  ! Passa alla prossima combinazione se la cartella non esiste
        end if

        ! Creazione del nome del file di output
        output_filename = trim(tracking_dir) // "SR1_" // trim(visco_short(v)) // "_" // &
                          trim(plast_short(p)) // "_" // trim(force_short(f)) // ".txt"

        ! Controllo se il file esiste già
        inquire(file=output_filename, exist=file_exists)
        if (file_exists) then
          print *, "File già esistente:", visco_short(v), plast_short(p), force_short(f)
          cycle  ! Salta alla prossima iterazione senza elaborare di nuovo
        else 
          print*, "Creo ",visco_short(v), plast_short(p), force_short(f)
        endif

        ! Apertura del file di output
        open(unit=10, file=output_filename, status='unknown')
        write(10,*) 'Time [Myr],  Nel,  Nsrmin,  SRmean, xmin, xmax'

        open(unit=30, file=trench_filename, status='old', iostat=stat)
          read(30, '(A)') line


        ! Loop sui file dei nodi per il tracking del trench
        do i = 0, 60  
          Nel=0 ! Numero di nodi all'interno della finestra
          Nemin = 0 ! Numero di nodi con sr >= srmin
          SRmean = 0.D0 ! Media dello strain rate per i nodi con sr >= srmin

          write(filename, '(A,I5.5,A)') 'Nodi.', i, '.txt'
          open(unit=20, file=trim(grid_dir)//trim(filename), status='old', iostat=stat)
          
          read(30,'(F10.3,2F12.3,3F12.3)') tT, xT
          if (stat.ne.0) then
            close(10)  
            exit  ! Esce se il file non esiste
          end if

          read(20, '(A)') line  

          ! Trova la prima cifra nella riga per estrarre il numero
          pos = scan(line, '0123456789')  
          if (pos > 0) then
            read(line(pos:), *) time
          else
            cycle
          end if

          do j = 1, 9  
            read(20,*)  
          end do


          do while (.true.)
            read(20,*, iostat=stat) node,x,y,a4,a5,a6,a7,a8,sr
            
            if (stat.ne.0) exit
            if (x.gt.xT-xmin.and.x.lt.xT+xmax.and.y.gt.ymin.and.y.lt.ymax) then
              
              Nel = Nel + 1
              
              ! Controllo per l'area deformata'
              if (sr >= srmin) then
                Nemin = Nemin + 1
                SRmean = SRmean + sr 
              end if
            end if

          end do
          close(20)

          write(10, '(F10.3, 2F12.3, 1E12.3, 2F12.3)') time, Nel, Nemin, SRmean/Nemin, xT - xmin, xT + xmax   

        end do

        close(10)
        close(30)
      end do
    end do
  end do

  print *, 'Tracking dello strain rate completato per tutti i modelli.'
end program sr_tracking
