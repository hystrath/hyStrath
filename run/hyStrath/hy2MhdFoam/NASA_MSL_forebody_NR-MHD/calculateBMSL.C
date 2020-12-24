#include <stdio.h>
#include <malloc.h>
#include <math.h>

void calculateCoilPoints(double** coil_points, const int steps, double x, double y, double z, double radius)
{
  //double center[3];
  double angle;
  const double pi = 3.1415926;
  int i;

  for (i = 0; i < steps; i++)
  {
    angle = 2*pi*i/(steps-1);
    //printf("%d %f\n", i, angle);
    coil_points[i][0] = x;
    coil_points[i][1] = y + radius*cos(angle);
    coil_points[i][2] = z + radius*sin(angle);
    //printf("%f %f %f \n", coil_points[i][0], coil_points[i][1], coil_points[i][2]);
    //fprintf(B, "%f %f %f \n", coil_points[0][i], coil_points[1][i], coil_points[2][i]);
  }
  return;
}
void allocate2dArray(double **a, int n, int m)
{
  int i, j;
  a = (double**)malloc(n*sizeof(double*));
  for (i = 0; i < n; i++)
     a[i] = (double*)malloc(m*sizeof(double));

  for (i = 0; i < n; i++)
  {
    for (j = 0; j < m; j++)
      a[i][j] = 0.0;
  }
}
void free2dArray(double **a, int n)
{
  int i;
  for(i = 0; i < n; i++)
  {
   free(a[i]);
  }
  free(a);
}


int main()
{
  FILE *ccx;
  FILE *ccy;
  FILE *ccz;
  int steps = 101;
  //double x = 0.1*(2/sqrt(3.0)-1);
  //double x = -0.38452994616;
  double x = 0.383191068;
  double y = 0.0;
  double z = 0.0;
  double radius = 0.31989224615;
  double** coil_points = NULL;
  double** cell_centers;
  FILE* B_file = fopen("0/B", "w");
  FILE* coil_file = fopen("coil", "w");
  FILE* centers_file = fopen("centers", "w");
  char buff[256];
  char buffx[256], buffy[256], buffz[256];
  int line_num_start = 21;
  int num_cells;
  int counter;

  double pi = 3.1415926;
  double mu0 = 4*pi*pow(10, -7);
  double B_target = 0.5*5.0/1.3; ///B at the center
  double I = 2*radius*B_target/mu0;
  double r[3], dl[3], r_unit[3];
  double r_mag;
  double c;
  double B[3];

  int i, j, k;

  //===allocate memory for coil points================
  coil_points = (double**)malloc(steps*sizeof(double*));
  for (i = 0; i < steps; i++)
     coil_points[i] = (double*)malloc(3*sizeof(double));

  for (i = 0; i < steps; i++)
  {
    for (j = 0; j < 3; j++)
      coil_points[i][j] = 0.0;
  }
  //===================================================
  calculateCoilPoints(coil_points, steps, x, y, z, radius);
  //===========read cell centers=======================
  ccx = fopen("0/Cx", "r");
  ccy = fopen("0/Cy", "r");
  ccz = fopen("0/Cz", "r");
  counter = 0;
  while (fgets(buff, sizeof(buff), ccx) != 0)
  {
    if (counter == line_num_start)
    {
      sscanf(buff, "%d", &num_cells);
      printf("%d\n", num_cells);
      break;
    }
    counter++;
  }

  //========allocate memory for cell centers array=========
  cell_centers = (double**)malloc(num_cells*sizeof(double*));
  for (i = 0; i < num_cells; i++)
     cell_centers[i] = (double*)malloc(3*sizeof(double));

  for (i = 0; i < num_cells; i++)
  {
    for (j = 0; j < 3; j++)
      cell_centers[i][j] = 0.0;
  }
  //=========================================================

  //=======Reading cell center coordinates into memory=======
  i = 0;
  while (fgets(buffx, sizeof(buff), ccx) != 0)
  {
    if (counter > line_num_start && counter < line_num_start+num_cells+1)
    {
      sscanf(buffx, "%lf", &cell_centers[i][0]);
      //printf("%f %f %f\n", cell_centers[i][0], cell_centers[i][1], cell_centers[i][2]);
      i++;
    }
    counter++;
  }


  i = 0;
  counter = 0;
  while (fgets(buffy, sizeof(buff), ccy) != 0)
  {
    if (counter > line_num_start+1 && counter < line_num_start+num_cells+2)
    {
      sscanf(buffy, "%lf", &cell_centers[i][1]);
      //printf("%f %f %f\n", cell_centers[i][0], cell_centers[i][1], cell_centers[i][2]);
      i++;
    }
    counter++;
  }

  i = 0;
  counter = 0;
  while (fgets(buffz, sizeof(buff), ccz) != 0)
  {
    if (counter > line_num_start+1 && counter < line_num_start+num_cells+2)
    {
      sscanf(buffz, "%lf", &cell_centers[i][2]);
      //printf("%f %f %f\n", cell_centers[i][0], cell_centers[i][1], cell_centers[i][2]);
      i++;
    }
    counter++;
  }
  //=========================================================

  //======Starting file...===================================
  fprintf(B_file, "/*--------------------------------*- C++ -*----------------------------------*\\\n\
| =========                 |                                                 |\n\
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n\
|  \\    /   O peration     | Version:  v1706                                 |\n\
|   \\\\  /    A nd           | Web:      www.OpenFOAM.com                      |\n\
|    \\\\/     M anipulation  |                                                 |\n\
\\*---------------------------------------------------------------------------*/\n\
FoamFile\n\
{\n\
    version     2.0;\n\
    format      ascii;\n\
    class       volVectorField;\n\
    location    \"0\";\n\
    object      B;\n\
}\n\
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n\
dimensions\t[1 0 -2 0 0 -1 0];\n\n\
internalField   nonuniform List<vector>\n");

  fprintf(B_file, "%d\n(\n", num_cells);

  //======Calculating B components===========================
  for (k = 0; k < num_cells; k++)
  {
    fprintf(centers_file, "%f %f %f;\n", cell_centers[k][0], cell_centers[k][1], cell_centers[k][2]);
    B[0] = 0.0;
    B[1] = 0.0;
    B[2] = 0.0;
    for (i = 0; i < steps-1; i++)
    {
      dl[0] = coil_points[i+1][0] - coil_points[i][0];
      dl[1] = coil_points[i+1][1] - coil_points[i][1];
      dl[2] = coil_points[i+1][2] - coil_points[i][2];
      r[0] = cell_centers[k][0] - coil_points[i][0];
      r[1] = cell_centers[k][1] - coil_points[i][1];
      r[2] = cell_centers[k][2] - coil_points[i][2];
      r_mag = sqrt(pow(r[0], 2) + pow(r[1], 2) + pow(r[2], 2));
      r_unit[0] = r[0]/r_mag;
      r_unit[1] = r[1]/r_mag;
      r_unit[2] = r[2]/r_mag;
      c = mu0*I/(4*pi*pow(r_mag, 2));
      B[0] += c*(dl[1]*r_unit[2] - dl[2]*r_unit[1]);
      B[1] += c*(dl[2]*r_unit[0] - dl[0]*r_unit[2]);
      B[2] += c*(dl[0]*r_unit[1] - dl[1]*r_unit[0]);

    }
    if (sqrt((pow(B[0], 2) + pow(B[1], 2) + pow(B[2], 2))) > 0.01)
      {
         printf("In cell #%d\n", k);
         printf("coordinates are: %f %f %f\n", cell_centers[k][0], cell_centers[k][1], cell_centers[k][2]);
         printf("B is: %d %f %f %f\n\n", k, B[0], B[1], B[2]);
      }
    //printf("%f\n", sqrt((pow(B[0], 2) + pow(B[1], 2) + pow(B[2], 2))));
    fprintf(B_file, "(%f %f %f)\n", B[0], B[1], B[2]); //writing B into file
    //printf("%d %f %f %f\n", k, B[0], B[1], B[2]);
  }
  printf("%f %f %f\n", cell_centers[0][0], cell_centers[0][1], cell_centers[0][2]);
  printf("%f %f %f\n", cell_centers[num_cells-1][0], cell_centers[num_cells-1][1], cell_centers[num_cells-1][2]);
  //===================================================

  //=====Finishing the file===========================

  fprintf(B_file, ")\n;\n\nboundaryField\n{\nfront\n{\n\ttype\t\twedge;\n}\nback\n{\n\ttype\t\twedge;\n}");
  fprintf(B_file, "\nobject\n{\n\ttype\t\tzeroGradient;\n}\n");
  fprintf(B_file, "inlet\n{\n\ttype\t\tzeroGradient;\n}\noutlet\n{\n\ttype\t\tzeroGradient;\n}\n}\n");
  fprintf(B_file, "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //");

  fprintf(coil_file, "[");
  for (i = 0; i < steps; i++)
  {
    fprintf(coil_file, "%f %f %f;\n", coil_points[i][0], coil_points[i][1], coil_points[i][2]);
  }
  fprintf(coil_file, "]");
  free2dArray(coil_points, steps);
  free2dArray(cell_centers, num_cells);
  fclose(B_file);
  fclose(coil_file);
  fclose(centers_file);
  fclose(ccx);
  fclose(ccy);
  fclose(ccz);
  return 0;
}


