# include "constant.c"

# define lx 2.
# define ly 1.
# define lz 1.
# define dx (lx/Nx)
# define dy (ly/Ny)
# define dz (lz/Nz)

# define norm(a,b,c) ((a)*(a)+(b)*(b)+(c)*(c))
# define modul(a,b,c) sqrt(norm(a,b,c))

//here is the main variables and structures of data
typedef double meshArr[Nx+2][Ny+2][Nz+2];


