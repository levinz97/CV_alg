#ifndef _LIBRARIES_H_
#define _LIBRARIES_H_

void read_pgm_header
(
    char  *filename,        /* in     : name of pgm file */
    long  *header_end_pos,  /* in+out : position of header end */
    int *nx,              /* out    : size in x-direction */
    int *ny               /* out    : size in y-direction */
)
/* reads PGM header */
{
    char row[80]; /* read buffer */
    FILE *file;   /* file pointer */

    /* try to open file */
    file = fopen(filename,"rb");

    /* if file not found */
    if (file==NULL)
    {
  printf("Reading %s FAILED",filename);
  printf("\n\n PROGRAM ABORTED !!! \n\n");
  exit(0);
    }

    /* read header */
    fgets (row, 300, file);
    fgets (row, 300, file);
    while (row[0]=='#') fgets(row, 300, file);

    /* read image dimensions */
    sscanf (row, "%d %d", nx, ny);
    fgets (row, 300, file);

    /* determine end of header */
    *header_end_pos=ftell(file);

    /* close file */
    fclose(file);
}

/* -------------------------------------------------------------------------- */

void read_pgm_data
(
    char   *filename,      /* in   : name of PGM file */
    long   data_start_pos, /* in   : position of data start */
    float **u,            /* out  : image */
    int  nx,             /* in   : size in x-direction */
    int  ny,             /* in   : size in y-direction */
    int  bx,             /* in   : boundary in x-direction */
    int  by              /* in   : boundary in y-direction */
)

/* reads PGM data */

{
    int i,j;  /* loop variables */
    FILE *file; /* file pointer */

    /* open file */
    file = fopen(filename,"rb");

    /* start at begin of data */
    fseek(file,data_start_pos,SEEK_SET);

    /* read image data */
    for (j=by; j<ny+by; j++)
  for (i=bx; i<nx+bx; i++)
  {
      u[i][j] = (float) getc (file);
  }

    /* close file */
    fclose(file);
}

/*--------------------------------------------------------------------------*/

void read_ppm_header
(
    char  *filename,        /* in     : name of pgm file */
    long  *header_end_pos,  /* in+out : position of header end */
    int *nx,              /* out    : size in x-direction */
    int *ny               /* out    : size in y-direction */
)
/* reads PGM header */
{
    char row[80]; /* read buffer */
    FILE *file;   /* file pointer */

    /* try to open file */
    file = fopen(filename,"rb");

    /* if file not found */
    if (file==NULL)
    {
  printf("Reading %s FAILED",filename);
  printf("\n\n PROGRAM ABORTED !!! \n\n");
  exit(0);
    }

    /* read header */
    fgets (row, 300, file);
    fgets (row, 300, file);
    while (row[0]=='#') fgets(row, 300, file);

    /* read image dimensions */
    sscanf (row, "%d %d", nx, ny);
    fgets (row, 300, file);

    /* determine end of header */
    *header_end_pos=ftell(file);

    /* close file */
    fclose(file);
}

/* -------------------------------------------------------------------------- */

void read_ppm_data_channelwise
(
    char   *filename,      /* in   : name of PGM file */
    long   data_start_pos, /* in   : position of data start */
    float **r,            /* out  : R channel of RGB image */
    float **g,            /* out  : G channel of RGB image */
    float **b,            /* out  : B channel of RGB image */
    int    nx,             /* in   : size in x-direction */
    int    ny,             /* in   : size in y-direction */
    int    bx,             /* in   : boundary in x-direction */
    int    by              /* in   : boundary in y-direction */
)

/* reads PGM data */

{
    int i,j;  /* loop variables */
    FILE *file; /* file pointer */

    /* open file */
    file = fopen(filename,"rb");

    /* start at begin of data */
    fseek(file,data_start_pos,SEEK_SET);

    /* read image data */
    for (j=by; j<ny+by; j++)
  for (i=bx; i<nx+bx; i++)
  {
      r[i][j] = (float) getc (file);
      g[i][j] = (float) getc (file);
      b[i][j] = (float) getc (file);
  }

    /* close file */
    fclose(file);
}

/* -------------------------------------------------------------------------- */

void write_pgm_data
(
    char *filename, /* in : file name */
    float **u,     /* in : image */
    int    nx,      /* in : size in x-direction */
    int    ny,      /* in : size in y-direction */
    int    bx,      /* in : boundary in x-direction */
    int    by       /* in : boundary in y-direction */
)

/* writes PGM data */

{

    FILE *file;  /* file pointer */
    int i,j;   /* loop variables */
    float help; /* tmp variable */
    unsigned char byte; /* variable for conversion */

    /* open file */
    file = fopen (filename, "ab");

    /* write image data */
    for (j=by; j<ny+by; j++)
  for (i=bx; i<nx+bx; i++)
  {
      help = u[i][j];
      if (help < 0.0)
    byte = (unsigned char)(0.0);
      else if (help > 255.0)
    byte = (unsigned char)(255.0);
      else
    byte = (unsigned char)(help);
      fwrite (&byte, sizeof(unsigned char), 1, file);
  }

    /* close file */
    fclose(file);
}

/* -------------------------------------------------------------------------- */
void write_ppm_data_channelwise
(
    char   *filename,      /* in   : name of PGM file */
    float **r,            /* out  : R channel of RGB image */
    float **g,            /* out  : G channel of RGB image */
    float **b,            /* out  : B channel of RGB image */
    int    nx,             /* in   : size in x-direction */
    int    ny,             /* in   : size in y-direction */
    int    bx,             /* in   : boundary in x-direction */
    int    by              /* in   : boundary in y-direction */
)

/* reads PPM data channelwise*/

{

    FILE *file;  /* file pointer */
    int   i,j;   /* loop variables */
    float help; /* tmp variable */
    unsigned char byte; /* variable for conversion */

    /* open file */
    file = fopen (filename, "ab");


    /* write image data */
    for (j=by; j<ny+by; j++)
  for (i=bx; i<nx+bx; i++)
  {
      help = r[i][j];
      if (help < 0.0)
    byte = (unsigned char)(0.0);
      else if (help > 255.0)
    byte = (unsigned char)(255.0);
      else
    byte = (unsigned char)(help);
      fwrite (&byte, sizeof(unsigned char), 1, file);

      help = g[i][j];
      if (help < 0.0)
    byte = (unsigned char)(0.0);
      else if (help > 255.0)
    byte = (unsigned char)(255.0);
      else
    byte = (unsigned char)(help);
      fwrite (&byte, sizeof(unsigned char), 1, file);

      help = b[i][j];
      if (help < 0.0)
    byte = (unsigned char)(0.0);
      else if (help > 255.0)
    byte = (unsigned char)(255.0);
      else
    byte = (unsigned char)(help);
      fwrite (&byte, sizeof(unsigned char), 1, file);
  }

    /* close file */
    fclose(file);
}

/*--------------------------------------------------------------------------*/

void alloc_vector

     (float **vector,   /* vector */
      long  n)          /* size */

     /* allocates storage for a vector of size n */


{
*vector = (float *) malloc (n * sizeof(float));
if (*vector == NULL)
   {
   printf("alloc_vector: not enough storage available\n");
   exit(1);
   }
return;
}

/*--------------------------------------------------------------------------*/

void alloc_matrix

     (float ***matrix,  /* matrix */
      long  nx,         /* size in x direction */
      long  ny)         /* size in y direction */

     /* allocates storage for matrix of size nx * ny */

{
long i;

*matrix = (float **) malloc (nx * sizeof(float *));
if (*matrix == NULL)
   {
   printf("alloc_matrix: not enough storage available\n");
   exit(1);
   }
for (i=0; i<nx; i++)
    {
    (*matrix)[i] = (float *) malloc (ny * sizeof(float));
    if ((*matrix)[i] == NULL)
       {
       printf("alloc_matrix: not enough storage available\n");
       exit(1);
       }
    }
return;
}

/*--------------------------------------------------------------------------*/

void disalloc_vector

     (float *vector,    /* vector */
      long  n)          /* size */

     /* disallocates storage for a vector of size n */

{
free(vector);
return;
}

/*--------------------------------------------------------------------------*/

void disalloc_matrix

     (float **matrix,   /* matrix */
      long  nx,         /* size in x direction */
      long  ny)         /* size in y direction */

     /* disallocates storage for matrix of size nx * ny */

{
long i;
for (i=0; i<nx; i++)
    free(matrix[i]);
free(matrix);
return;
}

/*--------------------------------------------------------------------------*/

void dummies

     (float **v,        /* image matrix */
      long  nx,         /* size in x direction */
      long  ny)         /* size in y direction */

/* creates dummy boundaries by mirroring */

{
long i, j;  /* loop variables */

for (i=1; i<=nx; i++)
    {
    v[i][0]    = v[i][1];
    v[i][ny+1] = v[i][ny];
    }

for (j=0; j<=ny+1; j++)
    {
    v[0][j]    = v[1][j];
    v[nx+1][j] = v[nx][j];
    }
return;
}

/*--------------------------------------------------------------------------*/

void presmooth

     (float    **f,       /* input: original image */
      float    **u,       /* output: smoothed */
      long     nx,        /* image dimension in x direction */
      long     ny,        /* image dimension in y direction */
      float    sigma)     /* standard deviation of Gaussian */



/*
 Gaussian convolution. Copyright by Joachim Weickert 5/2000
*/


{
long    i, j, p;              /* loop variables */
long    length;               /* convolution vector: 0..length */
float   sum;                  /* for summing up */
float   *conv;                /* convolution vector */
float   *help;                /* row or column with dummy boundaries */


/* ------------------------ diffusion in x direction -------------------- */

/* calculate length of convolution vector */
length = (long)(3.0 * sigma) + 1;
if (length > nx)
   {
   printf("gauss_conv: sigma too large \n");
   exit(0);
   }

/* allocate storage for convolution vector */
alloc_vector (&conv, length+1);

/* calculate entries of convolution vector */
for (i=0; i<=length; i++)
    conv[i] = 1 / (sigma * sqrt(2.0 * 3.1415926))
              * exp (- (i * i) / (2.0 * sigma * sigma));

/* normalization */
sum = conv[0];
for (i=1; i<=length; i++)
    sum = sum + 2.0 * conv[i];
for (i=0; i<=length; i++)
    conv[i] = conv[i] / sum;

/* allocate storage for a row */
alloc_vector (&help, nx+length+length);

for (j=1; j<=ny; j++)
    {
    /* copy in row vector */
    for (i=1; i<=nx; i++)
        help[i+length-1] = f[i][j];

    /* assign reflecting boundary conditions */
    for (p=1; p<=length; p++)
      {
  help[length-p]      = help[length+p-1];
  help[nx+length-1+p] = help[nx+length-p];
      }

    /* convolution step */
    for (i=length; i<=nx+length-1; i++)
        {
        /* calculate convolution */
        sum = conv[0] * help[i];
        for (p=1; p<=length; p++)
            sum = sum + conv[p] * (help[i+p] + help[i-p]);
        /* write back */
        u[i-length+1][j] = sum;
        }
    } /* for j */

/* disallocate storage for a row */
disalloc_vector (help, nx+length+length);

/* disallocate convolution vector */
disalloc_vector (conv, length + 1);


/* ------------------------ diffusion in y direction -------------------- */

/* calculate length of convolution vector */
length = (long)(3.0 * sigma) + 1;
if (length > ny)
   {
   printf("gauss_conv: sigma too large \n");
   exit(0);
   }

/* allocate storage for convolution vector */
alloc_vector (&conv, length + 1);

/* calculate entries of convolution vector */
for (j=0; j<=length; j++)
    conv[j] = 1 / (sigma * sqrt(2.0 * 3.1415927))
              * exp (- (j * j) / (2.0 * sigma * sigma));

/* normalization */
sum = conv[0];
for (j=1; j<=length; j++)
    sum = sum + 2.0 * conv[j];
for (j=0; j<=length; j++)
    conv[j] = conv[j] / sum;

/* allocate storage for a row */
alloc_vector (&help, ny+length+length);

for (i=1; i<=nx; i++)
    {
    /* copy in column vector */
    for (j=1; j<=ny; j++)
        help[j+length-1] = u[i][j];

    /* assign reflecting boundary conditions */
    for (p=1; p<=length; p++)
      {
  help[length-p]      = help[length+p-1];
  help[ny+length-1+p] = help[ny+length-p];
      }

    /* convolution step */
    for (j=length; j<=ny+length-1; j++)
        {
        /* calculate convolution */
        sum = conv[0] * help[j];
        for (p=1; p<=length; p++)
            sum = sum + conv[p] * (help[j+p] + help[j-p]);
        /* write back */
        u[i][j-length+1] = sum;
        }
    } /* for i */

/* disallocate storage for a row */
disalloc_vector (help, ny+length+length);

/* disallocate convolution vector */
disalloc_vector (conv, length+1);


return;

} /* gauss_conv */

#endif