  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!! FINITE ELEMENT METHOD !!!!!!!!!!!!!
  !!!!!!!!!!!!!     CST ELEMENTS      !!!!!!!!!!!!!
  !!!!!!!!!!!!!   MOHAMMAD SAJEDNIA   !!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program main
  implicit none
  integer :: nnd, nel, nbc, nload, ntrac, nvol, i, answer1, answer2, answer3
  integer, allocatable :: elements(:,:)
  real(8), dimension(:,:), allocatable :: nodes, nodal_loads, tractions, volume_loads, bc
  real :: E, nu, plane
  real(8) :: thickness
  real(8), allocatable :: D(:,:), K(:,:), F(:), U(:), reaction_forces(:)
  character(len=20) :: filename

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!             !!!!!!!!!!!!!
  !!!!!!!!!!!!! READ INPUTS !!!!!!!!!!!!!
  !!!!!!!!!!!!!             !!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  filename = 'nodes.txt'
  open(unit=10, file=filename, status='old', action='read')
  read(10,*) nnd
  allocate(nodes(nnd,2))
  do i=1,nnd
    read(10,*) nodes(i,:)
  end do
  close(10)

  filename = 'elements.txt'
  open(unit=11, file=filename, status='old', action='read')
  read(11,*) nel
  allocate(elements(nel,3))
  do i=1,nel
    read(11,*) elements(i,:)
  end do
  close(11)

  filename = 'boundaries.txt'
  open(unit=12, file=filename, status='old', action='read')
  read(12,*) nbc
  allocate(bc(nbc,3))
  do i=1,nbc
    read(12,*) bc(i,:)
  end do
  close(12)

  print*,'Is there any nodal loads? ( 1 for YES ) '
  read*, answer1
  if (answer1==1) then
  filename = 'nloads.txt'
  open(unit=13, file=filename, status='old', action='read')
  read(13,*) nload
  allocate(nodal_loads(nload,3))
  do i=1,nload
    read(13,*) nodal_loads(i,:)
  end do
  close(13)
  end if

  print*,'Is there any tractions? ( 1 for YES ) '
  read*, answer2
  if (answer2==1) then
  filename = 'tractions.txt'
  open(unit=14, file=filename, status='old', action='read')
  read(14,*) ntrac
  allocate(tractions(ntrac,4))
  do i=1,ntrac
    read(14,*) tractions(i,:)
  end do
  close(14)
  end if

  print*,'Is there any body forces? ( 1 For YES ) '
  read*, answer3
  if (answer3==1) then
  filename = 'vloads.txt'
  open(unit=15, file=filename, status='old', action='read')
  read(15,*) nvol
  allocate(volume_loads(nvol,3))
  do i=1,nvol
    read(15,*) volume_loads(i,:)
  end do
  close(15)
  end if

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!                       !!!!!!!!
  !!!!!!!! MATERIAL'S PROPERTIES !!!!!!!!
  !!!!!!!!                       !!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  E = 70e6
  nu = 0.25
  thickness = 0.02
  plane = 0 ! 0 if plane stress , 1 if plane strain

  if (plane==0) then
    allocate(D(3,3))
    D = (E/(1.0-nu**2.0)) * reshape([1.0, nu, 0.0, nu, 1.0, 0.0, 0.0, 0.0, (1.0-nu)/2.0], [3,3])
  elseif (plane==1) then
    allocate(D(3,3))
    D = (E/(1.0+nu)/(1.0-2.0*nu)) * reshape([1.0-nu, -nu, 0.0, -nu, 1.0-nu, 0.0, 0.0, 0.0, (1.0-2.0*nu)/2.0], [3,3])
  end if

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!                     !!!!!!!!
  !!!!!!!! CALLING SUBROUTINES !!!!!!!!
  !!!!!!!!                     !!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  allocate(F(2*nnd))
  F = 0.0
  call assemble_forces(thickness, nodes, answer1, nodal_loads, answer2, tractions, answer3, volume_loads, elements, F)

  allocate(K(2*nnd,2*nnd))
  K = 0.0
  call assemble_stiffness_matrix(nodes, nnd, elements, nel, thickness, D, K)

  allocate(U(2*nnd))
  U = 0.0
  call bcs_displacements(K, F, U, nbc, bc, nnd)

  allocate(reaction_forces(2*nnd))
  reaction_forces = 0.0
  call r_forces(K, U, bc, nbc, reaction_forces)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!                       !!!!!!!!
  !!!!!!!!    WRITING RESULTS    !!!!!!!!
  !!!!!!!!                       !!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  filename = 'Outputs.txt'
  open(unit=16, file=filename, status='replace', action='write')
  write(16,*) 'Nodal Displacements :'
  do i = 1, nnd
    write(16,*) 'Node', i, ': X =', U(2*i-1), ', Y =', U(2*i)
  end do
  write(16,*) 'Reaction Forces :'
  do i = 1, nnd
    write(16,*) 'Node', i, ': X =', reaction_forces(2*i-1), ', Y =', reaction_forces(2*i)
  end do
  close(16)

  print*, ''
  print*, 'Check Outputs.txt !'

  deallocate(nodes, elements, bc, D, K, F, U, reaction_forces)

  if (answer1==1) then
    deallocate(nodal_loads)
  end if

  if (answer2==1) then
    deallocate(tractions)
  end if

  if (answer3==1) then
    deallocate(volume_loads)
  end if

  stop

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!                          !!!!!!!!!!!!!
  !!!!!!!!!!!!! GAUSS ELIMINATION METHOD !!!!!!!!!!!!!
  !!!!!!!!!!!!!                          !!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine gauss(k,f,u)
  implicit none
  integer :: n
  real(8), dimension(:), intent(inout) :: u, f
  real(8), dimension(:,:), intent(inout) :: k
  real :: temp
  integer :: i, j, p, q
  real :: factor
  n = size(u)

  do i = 1, n-1
    p = i
    do j = i+1, n
      if (abs(k(j,i)) > abs(k(p,i))) then
        p = j
      end if
    end do

    if (p /= i) then
      do q = 1, n
        temp = k(i,q)
        k(i,q) = k(p,q)
        k(p,q) = temp
      end do
      temp = f(i)
      f(i) = f(p)
      f(p) = temp
    end if

    do j = i+1, n
      factor = k(j,i) / k(i,i)
      k(j,i) = 0.0
      do q = i+1, n
        k(j,q) = k(j,q) - factor*k(i,q)
      end do
      f(j) = f(j) - factor*f(i)
    end do
  end do

  u(n) = f(n) / k(n,n)
  do i = n-1, 1, -1
    u(i) = f(i)
    do j = i+1, n
      u(i) = u(i) - k(i,j)*u(j)
    end do
    u(i) = u(i) / k(i,i)
  end do

end subroutine gauss

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!                           !!!!!!!!!!!!!
  !!!!!!!!!!!!!  MATRICES MULTIPLICATION  !!!!!!!!!!!!!
  !!!!!!!!!!!!!                           !!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine sub_matmul(A, B, C, m, n, p)
  implicit none
  integer :: m, n, p
  real(8) :: A(m,n), B(n,p), C(m,p)
  integer :: i, j, k
  real :: sum0
  do i = 1, m
    do j = 1, p
      sum0 = 0.0
      do k = 1, n
        sum0 = sum0 + A(i,k) * B(k,j)
      C(i,j) = sum0
      end do
    end do
  end do
end subroutine sub_matmul

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!                          !!!!!!!!!!!!!
  !!!!!!!!!!!!!    MATRICES TRANSPOSE    !!!!!!!!!!!!!
  !!!!!!!!!!!!!                          !!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine sub_transpose(a, at)
  implicit none
  real(8), dimension(:,:), intent(in) :: a
  real(8), dimension(size(a, 2), size(a, 1)), intent(out) :: at
  integer :: i, j
  do i = 1, size(a, 1)
    do j = 1, size(a, 2)
      at(j, i) = a(i, j)
    end do
  end do
end subroutine sub_transpose

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!                     !!!!!!!!!!!!!
  !!!!!!!!!!!!!    FORCES VECTOR    !!!!!!!!!!!!!
  !!!!!!!!!!!!!                     !!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine assemble_forces(thickness, nodes, answer1, nodal_loads, answer2, tractions, answer3, vol_loads, elements, F)
    implicit none

    integer, intent(in) :: answer1, answer2, answer3
    real(8), intent(in) :: thickness
    real(8), dimension(:, :), intent(in) :: nodes
    real(8), dimension(:, :), intent(in) :: nodal_loads, tractions, vol_loads
    integer, dimension(:, :), intent(in) :: elements
    real(8), dimension(2*nnd, 1), intent(out) :: F

    integer :: i, j, ki, Fnode, first_node, end_node, vol_element, n1, n2, n3
    integer :: nload, ntrac, nvol,iii
    real(8) :: fx, fy, Dx, Dy, L, Tx, Ty, x1, x2, x3, y1, y2, y3, A, bx, by, pi
    real(8) :: Fmag, Fangle, traction_mag, traction_angle, vol_mag, vol_angle

    nload = size(nodal_loads,1)
    ntrac = size(tractions,1)
    nvol = size(vol_loads,1)
    pi = 3.14159265359

if (answer1 == 1) then
do i = 1, nload
    Fnode = nodal_loads(i, 1)
    Fmag = nodal_loads(i, 2)
    Fangle = nodal_loads(i, 3) / 180 * pi

    fx = Fmag * cos(Fangle)
    fy = Fmag * sin(Fangle)

    F(2 * Fnode - 1, 1) = F(2 * Fnode - 1, 1) + fx
    F(2 * Fnode, 1) = F(2 * Fnode, 1) + fy
end do
end if

if (answer2 == 1) then
do j = 1, ntrac
    first_node = tractions(j, 1)
    end_node = tractions(j, 2)
    traction_mag = tractions(j, 3)
    traction_angle = tractions(j, 4) / 180 * pi

    x1 = nodes(first_node, 1)
    x2 = nodes(end_node, 1)
    Dx = x2 - x1

    y1 = nodes(first_node, 2)
    y2 = nodes(end_node, 2)
    Dy = y2 - y1

    L = sqrt(Dx**2 + Dy**2)

    Tx = 0.5 * thickness * traction_mag * L * cos(traction_angle)
    Ty = 0.5 * thickness * traction_mag * L * sin(traction_angle)

    F(2 * first_node - 1, 1) = F(2 * first_node - 1, 1) + Tx
    F(2 * first_node, 1) = F(2 * first_node, 1) + Ty
    F(2 * end_node - 1, 1) = F(2 * end_node - 1, 1) + Tx
    F(2 * end_node, 1) = F(2 * end_node, 1) + Ty
end do
end if

if (answer3 == 1) then
do ki = 1, nvol
    vol_element = vol_loads(ki, 1)
    vol_mag = vol_loads(ki, 2)
    vol_angle = vol_loads(ki, 3) / 180 * pi

    n1 = elements(vol_element, 1)
    n2 = elements(vol_element, 2)
    n3 = elements(vol_element, 3)

    x1 = nodes(n1, 1)
    x2 = nodes(n2, 1)
    x3 = nodes(n3, 1)

    y1 = nodes(n1, 2)
    y2 = nodes(n2, 2)
    y3 = nodes(n3, 2)

    A = (x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2)) / 2.0

    bx = 1.0 / 3.0 * thickness * vol_mag * A * cos(vol_angle)
    by = 1.0 / 3.0 * thickness * vol_mag * A * sin(vol_angle)

    F(2 * n1 - 1, 1) = F(2 * n1 - 1, 1) + bx
    F(2 * n1, 1) = F(2 * n1, 1) + by
    F(2 * n2 - 1, 1) = F(2 * n2 - 1, 1) + bx
    F(2 * n2, 1) = F(2 * n2, 1) + by
    F(2 * n3 - 1, 1) = F(2 * n3 - 1, 1) + bx
    F(2 * n3, 1) = F(2 * n3, 1) + by
end do
end if

end subroutine assemble_forces

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!                               !!!!!!!!!!!!!
  !!!!!!!!!!!!!    GLOBAL STIFFNESS MATRIX    !!!!!!!!!!!!!
  !!!!!!!!!!!!!                               !!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine assemble_stiffness_matrix(nodes, nnd, elements, nel, thickness, D, K)

  implicit none

  integer, intent(in) :: nnd, nel
  real(8), dimension(:, :), intent(in) :: nodes
  integer, dimension(:, :), intent(in) :: elements
  real(8), intent(in) :: thickness
  real(8), dimension(3, 3), intent(in) :: D
  real(8), dimension(2*nnd, 2*nnd), intent(out) :: K

  integer :: i, n1, n2, n3, mm1, nn1, pp1, mm2, nn2, pp2
  real(8) :: x1, x2, x3, y1, y2, y3, A, m11, m21, m31, m12, m22, m32, m13, m23, m33
  real(8), dimension(3,6) :: B
  real(8), dimension(6,3) :: Bt, BTD
  real(8), dimension(6, 6) :: k1

  K = 0.0

  do i = 1, nel
    n1 = elements(i, 1)
    n2 = elements(i, 2)
    n3 = elements(i, 3)

    x1 = nodes(n1, 1)
    x2 = nodes(n2, 1)
    x3 = nodes(n3, 1)

    y1 = nodes(n1, 2)
    y2 = nodes(n2, 2)
    y3 = nodes(n3, 2)

    A = 0.5 * abs(x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2))

    m11 = (x2*y3 - x3*y2)/(2*A)
    m21 = (x3*y1 - x1*y3)/(2*A)
    m31 = (x1*y2 - y1*x2)/(2*A)
    m12 = (y2 - y3)/(2*A)
    m22 = (y3 - y1)/(2*A)
    m32 = (y1 - y2)/(2*A)
    m13 = (x3 - x2)/(2*A)
    m23 = (x1 - x3)/(2*A)
    m33 = (x2 - x1)/(2*A)

    B(1,1) = m12
    B(1,2) = 0
    B(1,3) = m22
    B(1,4) = 0
    B(1,5) = m32
    B(1,6) = 0
    B(2,1) = 0
    B(2,2) = m13
    B(2,3) = 0
    B(2,4) = m23
    B(2,5) = 0
    B(2,6) = m33
    B(3,1) = m13
    B(3,2) = m12
    B(3,3) = m23
    B(3,4) = m22
    B(3,5) = m33
    B(3,6) = m32

    call sub_transpose(B, Bt)

    mm1 = size(thickness * A * Bt,1)
    nn1 = size(thickness * A * Bt,2)
    pp1 = size(D,2)

    call sub_matmul(thickness * A * Bt, D, BTD, mm1, nn1, pp1)

    mm2 = size(BTD,1)
    nn2 = size(BTD,2)
    pp2 = size(B,2)

    call sub_matmul(BTD , B , k1, mm2, nn2, pp2)

    K(2*n1-1,2*n1-1) = K(2*n1-1,2*n1-1) + k1(1,1);
    K(2*n1-1,2*n1) = K(2*n1-1,2*n1) + k1(1,2);
    K(2*n1-1,2*n2-1) = K(2*n1-1,2*n2-1) + k1(1,3);
    K(2*n1-1,2*n2) = K(2*n1-1,2*n2) + k1(1,4);
    K(2*n1-1,2*n3-1) = K(2*n1-1,2*n3-1) + k1(1,5);
    K(2*n1-1,2*n3) = K(2*n1-1,2*n3) + k1(1,6);
    K(2*n1,2*n1-1) = K(2*n1,2*n1-1) + k1(2,1);
    K(2*n1,2*n1) = K(2*n1,2*n1) + k1(2,2);
    K(2*n1,2*n2-1) = K(2*n1,2*n2-1) + k1(2,3);
    K(2*n1,2*n2) = K(2*n1,2*n2) + k1(2,4);
    K(2*n1,2*n3-1) = K(2*n1,2*n3-1) + k1(2,5);
    K(2*n1,2*n3) = K(2*n1,2*n3) + k1(2,6);
    K(2*n2-1,2*n1-1) = K(2*n2-1,2*n1-1) + k1(3,1);
    K(2*n2-1,2*n1) = K(2*n2-1,2*n1) + k1(3,2);
    K(2*n2-1,2*n2-1) = K(2*n2-1,2*n2-1) + k1(3,3);
    K(2*n2-1,2*n2) = K(2*n2-1,2*n2) + k1(3,4);
    K(2*n2-1,2*n3-1) = K(2*n2-1,2*n3-1) + k1(3,5);
    K(2*n2-1,2*n3) = K(2*n2-1,2*n3) + k1(3,6);
    K(2*n2,2*n1-1) = K(2*n2,2*n1-1) + k1(4,1);
    K(2*n2,2*n1) = K(2*n2,2*n1) + k1(4,2);
    K(2*n2,2*n2-1) = K(2*n2,2*n2-1) + k1(4,3);
    K(2*n2,2*n2) = K(2*n2,2*n2) + k1(4,4);
    K(2*n2,2*n3-1) = K(2*n2,2*n3-1) + k1(4,5);
    K(2*n2,2*n3) = K(2*n2,2*n3) + k1(4,6);
    K(2*n3-1,2*n1-1) = K(2*n3-1,2*n1-1) + k1(5,1);
    K(2*n3-1,2*n1) = K(2*n3-1,2*n1) + k1(5,2);
    K(2*n3-1,2*n2-1) = K(2*n3-1,2*n2-1) + k1(5,3);
    K(2*n3-1,2*n2) = K(2*n3-1,2*n2) + k1(5,4);
    K(2*n3-1,2*n3-1) = K(2*n3-1,2*n3-1) + k1(5,5);
    K(2*n3-1,2*n3) = K(2*n3-1,2*n3) + k1(5,6);
    K(2*n3,2*n1-1) = K(2*n3,2*n1-1) + k1(6,1);
    K(2*n3,2*n1) = K(2*n3,2*n1) + k1(6,2);
    K(2*n3,2*n2-1) = K(2*n3,2*n2-1) + k1(6,3);
    K(2*n3,2*n2) = K(2*n3,2*n2) + k1(6,4);
    K(2*n3,2*n3-1) = K(2*n3,2*n3-1) + k1(6,5);
    K(2*n3,2*n3) = K(2*n3,2*n3) + k1(6,6);
  end do

end subroutine assemble_stiffness_matrix

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!                                             !!!!!!!!!!!!!
  !!!!!!!!!!!!!    BOUNDARY CONDITIONS AND DISPLACEMENTS    !!!!!!!!!!!!!
  !!!!!!!!!!!!!                                             !!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine bcs_displacements(K, F, U, nbc, bc, nnd)

  implicit none

  real(8), dimension(:,:), intent(in) :: K
  real(8), dimension(:), intent(inout) :: F
  real(8), dimension(:, :), intent(in) :: bc
  integer, intent(in) :: nnd, nbc
  real(8), dimension(:), intent(out) :: U

  integer :: i, j, node, direction, cnt
  real(8) :: value0
  integer, dimension(:), allocatable :: uu_zero
  real(8), dimension(:, :), allocatable :: K_F
  real(8), dimension(:), allocatable :: F_F, uf
  integer :: counter, cnti, cntj, cntii, adof

  adof = 2 * nnd - nbc ! Active degrees of freedom

  allocate(K_F(adof, adof))
  allocate(F_F(adof))
  allocate(uf(adof))
  allocate(uu_zero(nbc))

  U = 0.0d0
  cnt = 0

  do i = 1, nbc
    node = bc(i, 1)
    direction = bc(i, 2)
    value0 = bc(i, 3)

    U(2 * node - (2 - direction)) = value0

    F = F - K(:,2*node-(2-direction)) * value0

    cnt = cnt + 1
    uu_zero(cnt) = 2 * node - (2 - direction)

  end do

  K_F = 0.0d0
  cnti = 1

  do i = 1, 2*nnd
    if (all(uu_zero /= i)) then
      cntj = 1
      do j = 1, 2*nnd
        if (all(uu_zero /= j)) then
          K_F(cnti,cntj) = k(i, j)
          cntj = cntj + 1
        end if
      end do
      cnti = cnti + 1
    end if
  end do

  F_F = 0.0d0
  cntii = 1

  do i = 1, 2*nnd
    if (all(uu_zero /= i)) then
      F_F(cntii) = F(i)
      cntii = cntii + 1
    end if
  end do

  call gauss(K_F, F_F, uf)

  counter = 1
  do i = 1, 2 * nnd
    if (all(uu_zero /= i)) then
      U(i) = uf(counter)
      counter = counter + 1
    end if
  end do

  deallocate(K_F, F_F, uf, uu_zero)

end subroutine bcs_displacements

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!                       !!!!!!!!!!!!!
  !!!!!!!!!!!!!    REACTION FORCES    !!!!!!!!!!!!!
  !!!!!!!!!!!!!                       !!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine r_forces(K, U, bc, nbc, reaction_forces)
  implicit none
  real(8), intent(in) :: K(:,:)
  real(8), intent(in) :: U(:)
  real(8), intent(in) :: bc(:, :)
  integer, intent(in) :: nbc
  real(8), intent(out) :: reaction_forces(:)
  integer :: i, node0, direction0, mm, nn, pp

  reaction_forces = 0.0d0

  do i = 1, nbc
    node0 = bc(i, 1)
    direction0 = bc(i, 2)

  mm = 1
  nn = size(K,2)
  pp = 1

  call sub_matmul(K(2 * node0 - (2 - direction0), :), U, reaction_forces(2 * node0 - (2 - direction0):1), mm, nn, pp)

  end do

end subroutine r_forces

end program main

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!                       !!!!!!!!!!!!!
  !!!!!!!!!!!!!        THE END        !!!!!!!!!!!!!
  !!!!!!!!!!!!!                       !!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
