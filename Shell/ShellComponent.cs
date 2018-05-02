using System;
using System.Collections.Generic;

using Grasshopper.Kernel;
using Rhino.Geometry;

using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;

namespace Shell
{
    public class ShellComponent : GH_Component
    {
        public ShellComponent()
          : base("ShellCalculation", "SC",
              "Description",
              "Koala", "Shell")
        {
        }

        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddMeshParameter("IsoMesh", "IM", "The shell 6-node element isoparametric mesh, made by IsoMesher component", GH_ParamAccess.item);
            pManager.AddTextParameter("Boundary Conditions", "BDC", "Boundary Conditions in form x,y,z,vx,vy,vz,rx,ry,rz", GH_ParamAccess.list);
            pManager.AddTextParameter("Material properties", "Mat", "Material Properties", GH_ParamAccess.item, "210000,3600,4920000,4920000,79300,0.3,10");
            pManager.AddTextParameter("PointLoads", "PL", "Load given as Vector [N]", GH_ParamAccess.list);
            pManager.AddTextParameter("PointMoment", "PM", "Moment set in a point in [Nm]", GH_ParamAccess.list, "");
            pManager.AddBooleanParameter("Start calculations", "SC", "Set true to start calculations", GH_ParamAccess.item, false);
        }

        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddNumberParameter("Deformations", "Def", "Deformations", GH_ParamAccess.list);
            pManager.AddNumberParameter("Reactions", "R", "Reaction Forces", GH_ParamAccess.list);
            pManager.AddNumberParameter("Element stresses", "Strs", "The Stress in each element", GH_ParamAccess.list);
            pManager.AddNumberParameter("Element strains", "Strn", "The Strain in each element", GH_ParamAccess.list);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            #region Fetch inputs and assign to variables

            //Expected inputs
            Mesh mesh = new Mesh();                         //mesh in Mesh format
            List<MeshFace> faces = new List<MeshFace>();    //faces of mesh as a list
            List<Point3d> vertices = new List<Point3d>();   //vertices of mesh as a list

            List<string> bdctxt = new List<string>();       //Boundary conditions in string format
            List<string> loadtxt = new List<string>();      //loads in string format
            List<string> momenttxt = new List<string>();    //Moments in string format
            string mattxt = "";                             //Material in string format
            bool startCalc = false;

            if (!DA.GetData(0, ref mesh)) return;           //sets inputted mesh into variable
            if (!DA.GetDataList(1, bdctxt)) return;         //sets boundary conditions as string
            if (!DA.GetData(2, ref mattxt)) return;         //sets material properties as string
            if (!DA.GetDataList(3, loadtxt)) return;        //sets load as string
            if (!DA.GetDataList(4, momenttxt)) return;      //sets moment as string
            if (!DA.GetData(5, ref startCalc)) return;      //sets the boolean value for running the calculations

            if (!startCalc) return; //send return if startCalc is false

            foreach (var face in mesh.Faces)
            {
                faces.Add(face);
            }

            foreach (var vertice in mesh.Vertices)
            {
                vertices.Add(vertice);
            }
            #endregion

            //Interpret and set material parameters
            double E;       //Material Young's modulus, initial value 210000 [MPa]
            double A;       //Area for each element in same order as geometry, initial value CFS100x100 3600 [mm^2]
            double Iy;      //Moment of inertia about local y axis, initial value 4.92E6 [mm^4]
            double Iz;      //Moment of inertia about local z axis, initial value 4.92E6 [mm^4]
            double J;       //Polar moment of inertia
            double G;       //Shear modulus, initial value 79300 [mm^4]
            double nu;      //Poisson's ratio, initially 0.3
            double t;       //Thickness of shell
            SetMaterial(mattxt, out E, out A, out Iy, out Iz, out J, out G, out nu, out t);


            #region Prepares boundary conditions and loads for calculation

            //Interpret the BDC inputs (text) and create list of boundary condition (1/0 = free/clamped) for each dof.
            List<int> bdc_value = CreateBDCList(bdctxt, vertices);


            //Interpreting input load (text) and creating load list (do uble)
            List<double> load = CreateLoadList(loadtxt, momenttxt, vertices);
            #endregion

            #region Create global and reduced stiffness matrix
            //Create global stiffness matrix
            Matrix<double> K_tot = CreateGlobalStiffnessMatrix(faces, vertices, E, A, Iy, Iz, J, G, nu);

            //Create reduced K-matrix and reduced load list (removed free dofs)
            Matrix<double> K_red;
            Vector<double> load_red;
            CreateReducedGlobalStiffnessMatrix(bdc_value, K_tot, load, out K_red, out load_red);
            #endregion

        }

        private void CreateReducedGlobalStiffnessMatrix(List<int> bdc_value, Matrix<double> K, List<double> load, out Matrix<double> K_red, out Vector<double> load_red)
        {
            K_red = Matrix<double>.Build.DenseOfMatrix(K);
            List<double> load_redu = new List<double>(load);
            for (int i = 0, j = 0; i < load.Count; i++)
            {
                if (bdc_value[i] == 0)
                {
                    K_red = K_red.RemoveRow(i - j);
                    K_red = K_red.RemoveColumn(i - j);
                    load_redu.RemoveAt(i - j);
                    j++;
                }
            }
            load_red = Vector<double>.Build.DenseOfEnumerable(load_redu);
        }

        private static int GetGdofs(List<Point3d> vertices)
        {
            List<Point3d> uniqueNodes = new List<Point3d>();
            for (int i = 0; i < vertices.Count; i++)
            {
                Point3d tempNode = new Point3d(Math.Round(vertices[i].X, 2), Math.Round(vertices[i].Y, 2), Math.Round(vertices[i].Z, 2));
                if (!uniqueNodes.Contains(tempNode))
                {
                    uniqueNodes.Add(tempNode);
                }
            }
            return uniqueNodes.Count;
        }

        private Matrix<double> CreateGlobalStiffnessMatrix(List<MeshFace> faces, List<Point3d> vertices, double E, double A, double Iy, double Iz, double J, double G, double nu)
        {
            int ldof = 6;
            int gdofs = GetGdofs(vertices);
            var KG = Matrix<double>.Build.Dense(gdofs, gdofs);
            var KE = Matrix<double>.Build.Dense(3 * ldof, 3 * ldof);

            //temp t (until input is set up)
            double t = 100;

            for (int verticeIndex = 0; verticeIndex < vertices.Count / 3; verticeIndex += 3)
            {
                double y1 = vertices[verticeIndex].Y;
                double y2 = vertices[verticeIndex + 1].Y;
                double y3 = vertices[verticeIndex + 2].Y;

                double x1 = vertices[verticeIndex].X;
                double x2 = vertices[verticeIndex + 1].X;
                double x3 = vertices[verticeIndex + 2].X;

                double z1 = vertices[verticeIndex].Z;
                double z2 = vertices[verticeIndex + 1].Z;
                double z3 = vertices[verticeIndex + 2].Z;

                double[] xList = new double[3] { x1, x2, x3 };
                double[] yList = new double[3] { y1, y2, y3 };
                double[] zList = new double[3] { z1, z2, z3 };

                double area = 1 / 2 * Math.Sqrt(Math.Pow(x2 * y3 - x3 * y2, 2) + Math.Pow(x3 * y1 - x1 * y3, 2) + Math.Pow(x1 * y2 - x2 * y1, 2));

                Matrix<double> Ke = CreateElementStiffnessMatrix(xList, yList, zList, area, E, t, nu);
                



            int verticesPerFace = 3; //since triangle

            //Inputting values to correct entries in Global Stiffness Matrix
            for (int row = 0; row < KE.RowCount / verticesPerFace; row++)
            {
                for (int col = 0; col < KE.ColumnCount / verticesPerFace; col++)
                {
                    //top left 3x3 of K-element matrix
                    KG[verticeIndex * ldof + row, verticeIndex * ldof + col] += KE[row, col];
                    //top middle 3x3 of k-element matrix
                    KG[verticeIndex * ldof + row, (verticeIndex + 1) * ldof + col] += KE[row, col + ldof];
                    //top right 3x3 of k-element matrix  
                    KG[verticeIndex * ldof + row, (verticeIndex + 2) * ldof + col] += KE[row, col + ldof * 2];

                    //middle left 3x3 of k-element matrix
                    KG[(verticeIndex + 1) * ldof + row, verticeIndex * ldof + col] += KE[row + ldof, col];
                    //middle middle 3x3 of k-element matrix
                    KG[(verticeIndex + 1) * ldof + row, (verticeIndex + 1) * ldof + col] += KE[row + ldof, col + ldof];
                    //middle right 3x3 of k-element matrix
                    KG[(verticeIndex + 1) * ldof + row, (verticeIndex + 2) * ldof + col] += KE[row + ldof, col + ldof * 2];

                    //bottom left 3x3 of k-element matrix
                    KG[(verticeIndex + 2) * ldof + row, verticeIndex * ldof + col] += KE[row + ldof * 2, col];
                    //bottom middle 3x3 of k-element matrix
                    KG[(verticeIndex + 2) * ldof + row, (verticeIndex + 1) * ldof + col] += KE[row + ldof * 2, col + ldof];
                    //bottom right 3x3 of k-element matrix
                    KG[(verticeIndex + 2) * ldof + row, (verticeIndex + 2) * ldof + col] += KE[row + ldof * 2, col + ldof * 2];
                }
                //UNTESTED!!!
            }

            //Matrix<double> C = Matrix<double>.Build.DenseOfArray(new double[,] {
            //{  1/E, -nu/E, -nu/E},
            //{ },
            //{ },

            //    });
            //int eta = 0;
            //int eps = 0;

            #region old kmatrix code
            //int gdofs = points.Count * 6;
            //Matrix<double> K_tot = DenseMatrix.OfArray(new double[gdofs, gdofs]);

            //foreach (Line currentLine in geometry)
            //{
            //    double L = Math.Round(currentLine.Length, 6);

            //    Point3d p1 = new Point3d(Math.Round(currentLine.From.X, 2), Math.Round(currentLine.From.Y, 2), Math.Round(currentLine.From.Z, 2));
            //    Point3d p2 = new Point3d(Math.Round(currentLine.To.X, 2), Math.Round(currentLine.To.Y, 2), Math.Round(currentLine.To.Z, 2));

            //    double alpha = 0;

            //    double cx = (p2.X - p1.X) / L;
            //    double cy = (p2.Y - p1.Y) / L;
            //    double cz = (p2.Z - p1.Z) / L;
            //    double c1 = Math.Cos(alpha);
            //    double s1 = Math.Sin(alpha);
            //    double cxz = Math.Round(Math.Sqrt(Math.Pow(cx, 2) + Math.Pow(cz, 2)), 6);

            //    Matrix<double> gamma;

            //    if (Math.Round(cx, 6) == 0 && Math.Round(cz, 6) == 0)
            //    {
            //        gamma = Matrix<double>.Build.DenseOfArray(new double[,]
            //    {
            //        {      0, cy,  0},
            //        { -cy*c1,  0, s1},
            //        {  cy*s1,  0, c1},
            //    });
            //    }
            //    else
            //    {
            //        gamma = Matrix<double>.Build.DenseOfArray(new double[,]
            //    {
            //        {                     cx,       cy,                   cz},
            //        {(-cx*cy*c1 - cz*s1)/cxz,   cxz*c1,(-cy*cz*c1+cx*s1)/cxz},
            //        {   (cx*cy*s1-cz*c1)/cxz,  -cxz*s1, (cy*cz*s1+cx*c1)/cxz},
            //    });
            //    }

            //    var bd = Matrix<double>.Build;

            //    Matrix<double> T1;
            //    T1 = gamma.Append(bd.Dense(3, 9));
            //    Matrix<double> T2;
            //    T2 = bd.Dense(3, 3).Append(gamma);
            //    T2 = T2.Append(bd.Dense(3, 6));
            //    Matrix<double> T3;
            //    T3 = bd.Dense(3, 6).Append(gamma);
            //    T3 = T3.Append(bd.Dense(3, 3));
            //    Matrix<double> T4;
            //    T4 = bd.Dense(3, 9).Append(gamma);
            //    Matrix<double> T;
            //    T = T1.Stack(T2);
            //    T = T.Stack(T3);
            //    T = T.Stack(T4);

            //    //Matrix<double> T = SparseMatrix.OfArray(new double[,]
            //    //{
            //    //    { cx, cy, cz, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
            //    //    { cx, cy, cz, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
            //    //    { cx, cy, cz, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
            //    //    { 0, 0, 0, cx, cy, cz, 0, 0, 0, 0, 0, 0 },
            //    //    { 0, 0, 0, cx, cy, cz, 0, 0, 0, 0, 0, 0 },
            //    //    { 0, 0, 0, cx, cy, cz, 0, 0, 0, 0, 0, 0 },
            //    //    { 0, 0, 0, 0, 0, 0, cx, cy, cz, 0, 0, 0 },
            //    //    { 0, 0, 0, 0, 0, 0, cx, cy, cz, 0, 0, 0 },
            //    //    { 0, 0, 0, 0, 0, 0, cx, cy, cz, 0, 0, 0 },
            //    //    { 0, 0, 0, 0, 0, 0, 0, 0, 0, cx, cy, cz },
            //    //    { 0, 0, 0, 0, 0, 0, 0, 0, 0, cx, cy, cz },
            //    //    { 0, 0, 0, 0, 0, 0, 0, 0, 0, cx, cy, cz },
            //    //});

            //    Matrix<double> T_T = T.Transpose();

            //    double A1 = (E * A) / (L);

            //    double kz1 = (12 * E * Iz) / (L * L * L);
            //    double kz2 = (6 * E * Iz) / (L * L);
            //    double kz3 = (4 * E * Iz) / L;
            //    double kz4 = (2 * E * Iz) / L;

            //    double ky1 = (12 * E * Iy) / (L * L * L);
            //    double ky2 = (6 * E * Iy) / (L * L);
            //    double ky3 = (4 * E * Iy) / L;
            //    double ky4 = (2 * E * Iy) / L;

            //    double C1 = (G * J) / L;

            //    Matrix<double> K_elem = DenseMatrix.OfArray(new double[,]
            //    {
            //        { A1,    0,    0,    0,    0,    0,  -A1,    0,    0,    0,    0,    0 },
            //        {  0,  kz1,    0,    0,    0,  kz2,    0, -kz1,    0,    0,    0,  kz2 },
            //        {  0,    0,  ky1,    0, -ky2,    0,    0,    0, -ky1,    0, -ky2,    0 },
            //        {  0,    0,    0,   C1,    0,    0,    0,    0,    0,  -C1,    0,    0 },
            //        {  0,    0, -ky2,    0,  ky3,    0,    0,    0,  ky2,    0,  ky4,    0 },
            //        {  0,  kz2,    0,    0,    0,  kz3,    0, -kz2,    0,    0,    0,  kz4 },
            //        {-A1,    0,    0,    0,    0,    0,   A1,    0,    0,    0,    0,    0 },
            //        {  0, -kz1,    0,    0,    0, -kz2,    0,  kz1,    0,    0,    0, -kz2 },
            //        {  0,    0, -ky1,    0,  ky2,    0,    0,    0,  ky1,    0,  ky2,    0 },
            //        {  0,    0,    0,  -C1,    0,    0,    0,    0,    0,   C1,    0,    0 },
            //        {  0,    0, -ky2,    0,  ky4,    0,    0,    0,  ky2,    0,  ky3,    0 },
            //        {  0,  kz2,    0,    0,    0,  kz4,    0, -kz2,    0,    0,    0,  kz3 },
            //    });

            //    K_elem = K_elem.Multiply(T);
            //    K_elem = T_T.Multiply(K_elem);

            //    int node1 = points.IndexOf(p1);
            //    int node2 = points.IndexOf(p2);

            //    //System.Diagnostics.Debug.WriteLine("Node1: " + node1.ToString() + ", Node2: " + node2.ToString());

            //    //PrintMatrix(K_elem,"K_elem");

            //    //Inputting values to correct entries in Global Stiffness Matrix
            //    for (int i = 0; i < K_elem.RowCount / 2; i++)
            //    {
            //        //top left 3x3 of k-element matrix
            //        for (int j = 0; j < K_elem.ColumnCount / 2; j++)
            //        {
            //            K_tot[node1 * 6 + i, node1 * 6 + j] += K_elem[i, j];
            //        }
            //        //top right 3x3 of k-element matrix  
            //        for (int j = 0; j < K_elem.ColumnCount / 2; j++)
            //        {
            //            K_tot[node1 * 6 + i, node2 * 6 + j] += K_elem[i, j + 6];
            //        }
            //        //bottom left 3x3 of k-element matrix
            //        for (int j = 0; j < K_elem.ColumnCount / 2; j++)
            //        {
            //            K_tot[node2 * 6 + i, node1 * 6 + j] += K_elem[i + 6, j];
            //        }
            //        //bottom right 3x3 of k-element matrix
            //        for (int j = 0; j < K_elem.ColumnCount / 2; j++)
            //        {
            //            K_tot[node2 * 6 + i, node2 * 6 + j] += K_elem[i + 6, j + 6];
            //        }
            //    }
            //}
            #endregion

        }
            return KG;
            //NB! Consider calculating stresses and strains via this function to optimise calculation time!
            // el strain = Bq, el stress = DBq 
            //would be fastest to call calcstress&strain-method from this method since B is not saved outside this method!
        }

        private Matrix<double> CreateElementStiffnessMatrix(double[] xList, double[] yList, double[] zList, double Area, double E, double nu, double t)
        {
            //Må tranformeres til lokale koordinater!!
            double x1 = xList[0];
            double x2 = xList[1];
            double x3 = xList[2];

            double y1 = yList[0];
            double y2 = yList[1];
            double y3 = yList[2];

            double z1 = zList[0];
            double z2 = zList[1];
            double z3 = zList[2];

            double cosxX = -(x1 - x2) / Math.Pow((Math.Pow((x1 - x2), 2) + Math.Pow((y1 - y2), 2) + Math.Pow((z1 - z2), 2)), (1 / 2));
            double cosxY = -(y1 - y2) / Math.Pow((Math.Pow((x1 - x2), 2) + Math.Pow((y1 - y2), 2) + Math.Pow((z1 - z2), 2)), (1 / 2));
            double cosxZ = -(z1 - z2) / Math.Pow((Math.Pow((x1 - x2), 2) + Math.Pow((y1 - y2), 2) + Math.Pow((z1 - z2), 2)), (1 / 2));
            double cosyX = ((y1 - y2) * ((x1 - x2) * (y1 - y3) - (x1 - x3) * (y1 - y2)) + (z1 - z2) * ((x1 - x2) * (z1 - z3) - (x1 - x3) * (z1 - z2))) / Math.Pow((Math.Pow(((y1 - y2) * ((x1 - x2) * (y1 - y3) - (x1 - x3) * (y1 - y2)) + (z1 - z2) * ((x1 - x2) * (z1 - z3) - (x1 - x3) * (z1 - z2))), 2) + Math.Pow(((x1 - x2) * ((x1 - x2) * (y1 - y3) - (x1 - x3) * (y1 - y2)) - (z1 - z2) * ((y1 - y2) * (z1 - z3) - (y1 - y3) * (z1 - z2))), 2) + Math.Pow(((x1 - x2) * ((x1 - x2) * (z1 - z3) - (x1 - x3) * (z1 - z2)) + (y1 - y2) * ((y1 - y2) * (z1 - z3) - (y1 - y3) * (z1 - z2))), 2)), (1 / 2));
            double cosyY = -((x1 - x2) * ((x1 - x2) * (y1 - y3) - (x1 - x3) * (y1 - y2)) - (z1 - z2) * ((y1 - y2) * (z1 - z3) - (y1 - y3) * (z1 - z2))) / Math.Pow((Math.Pow(((y1 - y2) * ((x1 - x2) * (y1 - y3) - (x1 - x3) * (y1 - y2)) + (z1 - z2) * ((x1 - x2) * (z1 - z3) - (x1 - x3) * (z1 - z2))), 2) + Math.Pow(((x1 - x2) * ((x1 - x2) * (y1 - y3) - (x1 - x3) * (y1 - y2)) - (z1 - z2) * ((y1 - y2) * (z1 - z3) - (y1 - y3) * (z1 - z2))), 2) + Math.Pow(((x1 - x2) * ((x1 - x2) * (z1 - z3) - (x1 - x3) * (z1 - z2)) + (y1 - y2) * ((y1 - y2) * (z1 - z3) - (y1 - y3) * (z1 - z2))), 2)), (1 / 2));
            double cosyZ = -((x1 - x2) * ((x1 - x2) * (z1 - z3) - (x1 - x3) * (z1 - z2)) + (y1 - y2) * ((y1 - y2) * (z1 - z3) - (y1 - y3) * (z1 - z2))) / Math.Pow((Math.Pow(((y1 - y2) * ((x1 - x2) * (y1 - y3) - (x1 - x3) * (y1 - y2)) + (z1 - z2) * ((x1 - x2) * (z1 - z3) - (x1 - x3) * (z1 - z2))), 2) + Math.Pow(((x1 - x2) * ((x1 - x2) * (y1 - y3) - (x1 - x3) * (y1 - y2)) - (z1 - z2) * ((y1 - y2) * (z1 - z3) - (y1 - y3) * (z1 - z2))), 2) + Math.Pow(((x1 - x2) * ((x1 - x2) * (z1 - z3) - (x1 - x3) * (z1 - z2)) + (y1 - y2) * ((y1 - y2) * (z1 - z3) - (y1 - y3) * (z1 - z2))), 2)), (1 / 2));
            double coszX = ((y1 - y2) * (z1 - z3) - (y1 - y3) * (z1 - z2)) / Math.Pow((Math.Pow(((x1 - x2) * (y1 - y3) - (x1 - x3) * (y1 - y2)), 2) + Math.Pow(((x1 - x2) * (z1 - z3) - (x1 - x3) * (z1 - z2)), 2) + Math.Pow(((y1 - y2) * (z1 - z3) - (y1 - y3) * (z1 - z2)), 2)), (1 / 2));
            double coszY = -((x1 - x2) * (z1 - z3) - (x1 - x3) * (z1 - z2)) / Math.Pow((Math.Pow(((x1 - x2) * (y1 - y3) - (x1 - x3) * (y1 - y2)), 2) + Math.Pow(((x1 - x2) * (z1 - z3) - (x1 - x3) * (z1 - z2)), 2) + Math.Pow(((y1 - y2) * (z1 - z3) - (y1 - y3) * (z1 - z2)), 2)), (1 / 2));
            double coszZ = ((x1 - x2) * (y1 - y3) - (x1 - x3) * (y1 - y2)) / Math.Pow((Math.Pow(((x1 - x2) * (y1 - y3) - (x1 - x3) * (y1 - y2)), 2) + Math.Pow(((x1 - x2) * (z1 - z3) - (x1 - x3) * (z1 - z2)), 2) + Math.Pow(((y1 - y2) * (z1 - z3) - (y1 - y3) * (z1 - z2)), 2)), (1 / 2));

            Matrix<double> tf = Matrix<double>.Build.Dense(3, 3);
            tf[0, 0] = cosxX;
            tf[0, 1] = cosxY;
            tf[0, 2] = cosxZ;
            tf[1, 0] = cosyX;
            tf[1, 1] = cosyY;
            tf[1, 2] = cosyZ;
            tf[2, 0] = coszX;
            tf[2, 1] = coszY;
            tf[2, 2] = coszZ;

            Matrix<double> one = Matrix<double>.Build.Dense(1,1,1);
            tf = tf.DiagonalStack(one);
            var T = tf.DiagonalStack(tf);
            T = T.DiagonalStack(tf);
            Matrix<double> T_T = T.Transpose();

            Matrix<double> lcoord = Matrix<double>.Build.DenseOfArray(new double[,]
            {
                { x1, x2, x3 },
                { y1, y2, y3 },
                { z1, z2, z3 }
            });

            lcoord = T.Multiply(lcoord);
            x1 = lcoord[0, 0];
            x2 = lcoord[0, 1];
            x3 = lcoord[0, 2];
            y1 = lcoord[1, 0];
            y1 = lcoord[1, 1];
            y1 = lcoord[1, 2];

            double x13 = x1 - x3;
            double x32 = x3 - x2;
            double y23 = y2 - y3;
            double y31 = y3 - y1;

            double[] ga = new double[3];
            double[] my = new double[3];
            double[] a = new double[3];

            for (int i = 0; i < 4; i++)
            {
                double c, s;
                double len = Math.Sqrt(Math.Pow(xList[i + 1] - xList[i], 2) + Math.Pow(yList[i + 1] - yList[i], 2));
                if (xList[i+1]>xList[i])
                {
                    c = (xList[i + 1] - xList[i]) / len;
                    s = (yList[i + 1] - yList[i]) / len;
                }
                else if (xList[i+1]<xList[i])
                {
                    c = (xList[i] - xList[i + 1]) / len;
                    s = (yList[i] - yList[i + 1]) / len;
                }
                else
                {
                    c = 0;
                    s = 1;
                }
                ga[i] = (c * x32 - s * y23) / (2 * Area);
                my[i] = (c * x13 - s * y31) / (2 * Area);
                a[i] = ga[i] + my[i];
            }

            double ga4 = ga[0];
            double ga5 = ga[1];
            double ga6 = ga[2];
            double my4 = my[0];
            double my5 = my[1];
            double my6 = my[2];
            double a4 = a[0];
            double a5 = a[1];
            double a6 = a[2];

            Matrix<double> ke = Matrix<double>.Build.Dense(12, 12);
            //ke calculated in matlab script Simplest_shell_triangle.m in local xy coordinates
            ke[0, 0] = -(Area * E * t * (Math.Pow(x2, 2) - 4 * y2 * y3 - nu * Math.Pow(x2, 2) - nu * Math.Pow(x3, 2) - 2 * x2 * x3 + Math.Pow(x3, 2) + 2 * Math.Pow(y2, 2) + 2 * Math.Pow(y3, 2) + 2 * nu * x2 * x3)) / (2 * (Math.Pow(nu, 2) - 1) * Math.Pow((x1 * y2 - x2 * y1 - x1 * y3 + x3 * y1 + x2 * y3 - x3 * y2), 2));
            ke[0, 1] = (Area * E * t * (x2 - x3) * (y2 - y3) * (nu + 1)) / (2 * (Math.Pow(nu, 2) - 1) * Math.Pow((x1 * y2 - x2 * y1 - x1 * y3 + x3 * y1 + x2 * y3 - x3 * y2), 2));
            ke[0, 4] = (Area * E * t * (x1 * x2 - x1 * x3 - x2 * x3 + 2 * y1 * y2 - 2 * y1 * y3 - 2 * y2 * y3 - nu * Math.Pow(x3, 2) + Math.Pow(x3, 2) + 2 * Math.Pow(y3, 2) - nu * x1 * x2 + nu * x1 * x3 + nu * x2 * x3)) / (2 * (Math.Pow(nu, 2) - 1) * Math.Pow((x1 * y2 - x2 * y1 - x1 * y3 + x3 * y1 + x2 * y3 - x3 * y2), 2));
            ke[0, 5] = -(Area * E * t * (x2 * y1 - x3 * y1 - x2 * y3 + x3 * y3 + 2 * nu * x1 * y2 - nu * x2 * y1 - 2 * nu * x1 * y3 + nu * x3 * y1 + nu * x2 * y3 - 2 * nu * x3 * y2 + nu * x3 * y3)) / (2 * (Math.Pow(nu, 2) - 1) * Math.Pow((x1 * y2 - x2 * y1 - x1 * y3 + x3 * y1 + x2 * y3 - x3 * y2), 2));
            ke[0, 8] = -(Area * E * t * (x1 * x2 - x1 * x3 + x2 * x3 + 2 * y1 * y2 - 2 * y1 * y3 + 2 * y2 * y3 + nu * Math.Pow(x2, 2) - Math.Pow(x2, 2) - 2 * Math.Pow(y2, 2) - nu * x1 * x2 + nu * x1 * x3 - nu * x2 * x3)) / (2 * (Math.Pow(nu, 2) - 1) * Math.Pow((x1 * y2 - x2 * y1 - x1 * y3 + x3 * y1 + x2 * y3 - x3 * y2), 2));
            ke[0, 9] = -(Area * E * t * (x2 * y2 - x2 * y1 + x3 * y1 - x3 * y2 - 2 * nu * x1 * y2 + nu * x2 * y1 + 2 * nu * x1 * y3 + nu * x2 * y2 - nu * x3 * y1 - 2 * nu * x2 * y3 + nu * x3 * y2)) / (2 * (Math.Pow(nu, 2) - 1) * Math.Pow((x1 * y2 - x2 * y1 - x1 * y3 + x3 * y1 + x2 * y3 - x3 * y2), 2));
            ke[1, 0] = (Area * E * t * (x2 - x3) * (y2 - y3) * (nu + 1)) / (2 * (Math.Pow(nu, 2) - 1) * Math.Pow((x1 * y2 - x2 * y1 - x1 * y3 + x3 * y1 + x2 * y3 - x3 * y2), 2));
            ke[1, 1] = -(Area * E * t * (2 * Math.Pow(x2, 2) - 2 * y2 * y3 - nu * Math.Pow(y2, 2) - nu * Math.Pow(y3, 2) - 4 * x2 * x3 + 2 * Math.Pow(x3, 2) + Math.Pow(y2, 2) + Math.Pow(y3, 2) + 2 * nu * y2 * y3)) / (2 * (Math.Pow(nu, 2) - 1) * Math.Pow((x1 * y2 - x2 * y1 - x1 * y3 + x3 * y1 + x2 * y3 - x3 * y2), 2));
            ke[1, 4] = -(Area * E * t * (x1 * y2 - x1 * y3 - x3 * y2 + x3 * y3 - nu * x1 * y2 + 2 * nu * x2 * y1 + nu * x1 * y3 - 2 * nu * x3 * y1 - 2 * nu * x2 * y3 + nu * x3 * y2 + nu * x3 * y3)) / (2 * (Math.Pow(nu, 2) - 1) * Math.Pow((x1 * y2 - x2 * y1 - x1 * y3 + x3 * y1 + x2 * y3 - x3 * y2), 2));
            ke[1, 5] = (Area * E * t * (2 * x1 * x2 - 2 * x1 * x3 - 2 * x2 * x3 + y1 * y2 - y1 * y3 - y2 * y3 - nu * Math.Pow(y3, 2) + 2 * Math.Pow(x3, 2) + Math.Pow(y3, 2) - nu * y1 * y2 + nu * y1 * y3 + nu * y2 * y3)) / (2 * (Math.Pow(nu, 2) - 1) * Math.Pow((x1 * y2 - x2 * y1 - x1 * y3 + x3 * y1 + x2 * y3 - x3 * y2), 2));
            ke[1, 8] = -(Area * E * t * (x1 * y3 - x1 * y2 + x2 * y2 - x2 * y3 + nu * x1 * y2 - 2 * nu * x2 * y1 - nu * x1 * y3 + nu * x2 * y2 + 2 * nu * x3 * y1 + nu * x2 * y3 - 2 * nu * x3 * y2)) / (2 * (Math.Pow(nu, 2) - 1) * Math.Pow((x1 * y2 - x2 * y1 - x1 * y3 + x3 * y1 + x2 * y3 - x3 * y2), 2));
            ke[1, 9] = -(Area * E * t * (2 * x1 * x2 - 2 * x1 * x3 + 2 * x2 * x3 + y1 * y2 - y1 * y3 + y2 * y3 + nu * Math.Pow(y2, 2) - 2 * Math.Pow(x2, 2) - Math.Pow(y2, 2) - nu * y1 * y2 + nu * y1 * y3 - nu * y2 * y3)) / (2 * (Math.Pow(nu, 2) - 1) * Math.Pow((x1 * y2 - x2 * y1 - x1 * y3 + x3 * y1 + x2 * y3 - x3 * y2), 2));
            ke[2, 2] = -(E * Math.Pow(t, 3) * ((2 * Math.Pow(x32, 2) - (2 * ga6 * x13 * x32) / my6) * (2 * Math.Pow(x32, 2) + nu * (2 * Math.Pow(y23, 2) - (2 * ga6 * y23 * y31) / my6) - (2 * ga6 * x13 * x32) / my6) - (nu / 2 - 1 / 2) * Math.Pow((4 * x32 * y23 - (ga6 * (2 * x13 * y23 + 2 * x32 * y31)) / my6), 2) + (2 * Math.Pow(y23, 2) - (2 * ga6 * y23 * y31) / my6) * (2 * Math.Pow(y23, 2) + nu * (2 * Math.Pow(x32, 2) - (2 * ga6 * x13 * x32) / my6) - (2 * ga6 * y23 * y31) / my6))) / (192 * Math.Pow(Area, 3) * (Math.Pow(nu, 2) - 1));
            ke[2, 6] = -(E * Math.Pow(t, 3) * ((2 * Math.Pow(x13, 2) - (2 * my5 * x13 * x32) / ga5) * (2 * Math.Pow(x32, 2) + nu * (2 * Math.Pow(y23, 2) - (2 * ga6 * y23 * y31) / my6) - (2 * ga6 * x13 * x32) / my6) + (2 * Math.Pow(y31, 2) - (2 * my5 * y23 * y31) / ga5) * (2 * Math.Pow(y23, 2) + nu * (2 * Math.Pow(x32, 2) - (2 * ga6 * x13 * x32) / my6) - (2 * ga6 * y23 * y31) / my6) - (nu / 2 - 1 / 2) * (4 * x13 * y31 - (my5 * (2 * x13 * y23 + 2 * x32 * y31)) / ga5) * (4 * x32 * y23 - (ga6 * (2 * x13 * y23 + 2 * x32 * y31)) / my6))) / (192 * Math.Pow(Area, 3) * (Math.Pow(nu, 2) - 1));
            ke[2, 7] = -(E * Math.Pow(t, 3) * ((2 * x13 * x32 * (2 * Math.Pow(x32, 2) + nu * (2 * Math.Pow(y23, 2) - (2 * ga6 * y23 * y31) / my6) - (2 * ga6 * x13 * x32) / my6)) / ga5 + (2 * y23 * y31 * (2 * Math.Pow(y23, 2) + nu * (2 * Math.Pow(x32, 2) - (2 * ga6 * x13 * x32) / my6) - (2 * ga6 * y23 * y31) / my6)) / ga5 - ((2 * x13 * y23 + 2 * x32 * y31) * (nu / 2 - 1 / 2) * (4 * x32 * y23 - (ga6 * (2 * x13 * y23 + 2 * x32 * y31)) / my6)) / ga5)) / (192 * Math.Pow(Area, 3) * (Math.Pow(nu, 2) - 1));
            ke[2, 10] = -(E * Math.Pow(t, 3) * (2 * x13 * x32 * (a5 / ga5 + a6 / my6) * (2 * Math.Pow(x32, 2) + nu * (2 * Math.Pow(y23, 2) - (2 * ga6 * y23 * y31) / my6) - (2 * ga6 * x13 * x32) / my6) - (2 * x13 * y23 + 2 * x32 * y31) * (nu / 2 - 1 / 2) * (4 * x32 * y23 - (ga6 * (2 * x13 * y23 + 2 * x32 * y31)) / my6) * (a5 / ga5 + a6 / my6) + 2 * y23 * y31 * (a5 / ga5 + a6 / my6) * (2 * Math.Pow(y23, 2) + nu * (2 * Math.Pow(x32, 2) - (2 * ga6 * x13 * x32) / my6) - (2 * ga6 * y23 * y31) / my6))) / (192 * Math.Pow(Area, 3) * (Math.Pow(nu, 2) - 1));
            ke[2, 11] = -(E * Math.Pow(t, 3) * ((4 * x13 * x32 * (my6 * Math.Pow(x32, 2) + my6 * nu * Math.Pow(y23, 2) - ga6 * x13 * x32 - ga6 * nu * y23 * y31)) / Math.Pow(my6, 2) + (4 * y23 * y31 * (my6 * Math.Pow(y23, 2) + my6 * nu * Math.Pow(x32, 2) - ga6 * y23 * y31 - ga6 * nu * x13 * x32)) / Math.Pow(my6, 2) + (2 * (2 * x13 * y23 + 2 * x32 * y31) * (nu / 2 - 1 / 2) * (ga6 * x13 * y23 + ga6 * x32 * y31 - 2 * my6 * x32 * y23)) / Math.Pow(my6, 2))) / (192 * Math.Pow(Area, 3) * (Math.Pow(nu, 2) - 1));
            ke[4, 0] = (Area * E * t * (x1 * x2 - x1 * x3 - x2 * x3 + 2 * y1 * y2 - 2 * y1 * y3 - 2 * y2 * y3 - nu * Math.Pow(x3, 2) + Math.Pow(x3, 2) + 2 * Math.Pow(y3, 2) - nu * x1 * x2 + nu * x1 * x3 + nu * x2 * x3)) / (2 * (Math.Pow(nu, 2) - 1) * Math.Pow((x1 * y2 - x2 * y1 - x1 * y3 + x3 * y1 + x2 * y3 - x3 * y2), 2));
            ke[4, 1] = -(Area * E * t * (x1 * y2 - x1 * y3 - x3 * y2 + x3 * y3 - nu * x1 * y2 + 2 * nu * x2 * y1 + nu * x1 * y3 - 2 * nu * x3 * y1 - 2 * nu * x2 * y3 + nu * x3 * y2 + nu * x3 * y3)) / (2 * (Math.Pow(nu, 2) - 1) * Math.Pow((x1 * y2 - x2 * y1 - x1 * y3 + x3 * y1 + x2 * y3 - x3 * y2), 2));
            ke[4, 4] = -(Area * E * t * (Math.Pow(x1, 2) - 4 * y1 * y3 - nu * Math.Pow(x1, 2) - nu * Math.Pow(x3, 2) - 2 * x1 * x3 + Math.Pow(x3, 2) + 2 * Math.Pow(y1, 2) + 2 * Math.Pow(y3, 2) + 2 * nu * x1 * x3)) / (2 * (Math.Pow(nu, 2) - 1) * Math.Pow((x1 * y2 - x2 * y1 - x1 * y3 + x3 * y1 + x2 * y3 - x3 * y2), 2));
            ke[4, 5] = (Area * E * t * (x1 - x3) * (y1 - y3) * (nu + 1)) / (2 * (Math.Pow(nu, 2) - 1) * Math.Pow((x1 * y2 - x2 * y1 - x1 * y3 + x3 * y1 + x2 * y3 - x3 * y2), 2));
            ke[4, 8] = -(Area * E * t * (x1 * x2 + x1 * x3 - x2 * x3 + 2 * y1 * y2 + 2 * y1 * y3 - 2 * y2 * y3 + nu * Math.Pow(x1, 2) - Math.Pow(x1, 2) - 2 * Math.Pow(y1, 2) - nu * x1 * x2 - nu * x1 * x3 + nu * x2 * x3)) / (2 * (Math.Pow(nu, 2) - 1) * Math.Pow((x1 * y2 - x2 * y1 - x1 * y3 + x3 * y1 + x2 * y3 - x3 * y2), 2));
            ke[4, 9] = -(Area * E * t * (x1 * y1 - x1 * y2 - x3 * y1 + x3 * y2 + nu * x1 * y1 + nu * x1 * y2 - 2 * nu * x2 * y1 - 2 * nu * x1 * y3 + nu * x3 * y1 + 2 * nu * x2 * y3 - nu * x3 * y2)) / (2 * (Math.Pow(nu, 2) - 1) * Math.Pow((x1 * y2 - x2 * y1 - x1 * y3 + x3 * y1 + x2 * y3 - x3 * y2), 2));
            ke[5, 0] = -(Area * E * t * (x2 * y1 - x3 * y1 - x2 * y3 + x3 * y3 + 2 * nu * x1 * y2 - nu * x2 * y1 - 2 * nu * x1 * y3 + nu * x3 * y1 + nu * x2 * y3 - 2 * nu * x3 * y2 + nu * x3 * y3)) / (2 * (Math.Pow(nu, 2) - 1) * Math.Pow((x1 * y2 - x2 * y1 - x1 * y3 + x3 * y1 + x2 * y3 - x3 * y2), 2));
            ke[5, 1] = (Area * E * t * (2 * x1 * x2 - 2 * x1 * x3 - 2 * x2 * x3 + y1 * y2 - y1 * y3 - y2 * y3 - nu * Math.Pow(y3, 2) + 2 * Math.Pow(x3, 2) + Math.Pow(y3, 2) - nu * y1 * y2 + nu * y1 * y3 + nu * y2 * y3)) / (2 * (Math.Pow(nu, 2) - 1) * Math.Pow((x1 * y2 - x2 * y1 - x1 * y3 + x3 * y1 + x2 * y3 - x3 * y2), 2));
            ke[5, 4] = (Area * E * t * (x1 - x3) * (y1 - y3) * (nu + 1)) / (2 * (Math.Pow(nu, 2) - 1) * Math.Pow((x1 * y2 - x2 * y1 - x1 * y3 + x3 * y1 + x2 * y3 - x3 * y2), 2));
            ke[5, 5] = -(Area * E * t * (2 * Math.Pow(x1, 2) - 2 * y1 * y3 - nu * Math.Pow(y1, 2) - nu * Math.Pow(y3, 2) - 4 * x1 * x3 + 2 * Math.Pow(x3, 2) + Math.Pow(y1, 2) + Math.Pow(y3, 2) + 2 * nu * y1 * y3)) / (2 * (Math.Pow(nu, 2) - 1) * Math.Pow((x1 * y2 - x2 * y1 - x1 * y3 + x3 * y1 + x2 * y3 - x3 * y2), 2));
            ke[5, 8] = -(Area * E * t * (x1 * y1 - x2 * y1 - x1 * y3 + x2 * y3 + nu * x1 * y1 - 2 * nu * x1 * y2 + nu * x2 * y1 + nu * x1 * y3 - 2 * nu * x3 * y1 - nu * x2 * y3 + 2 * nu * x3 * y2)) / (2 * (Math.Pow(nu, 2) - 1) * Math.Pow((x1 * y2 - x2 * y1 - x1 * y3 + x3 * y1 + x2 * y3 - x3 * y2), 2));
            ke[5, 9] = -(Area * E * t * (2 * x1 * x2 + 2 * x1 * x3 - 2 * x2 * x3 + y1 * y2 + y1 * y3 - y2 * y3 + nu * Math.Pow(y1, 2) - 2 * Math.Pow(x1, 2) - Math.Pow(y1, 2) - nu * y1 * y2 - nu * y1 * y3 + nu * y2 * y3)) / (2 * (Math.Pow(nu, 2) - 1) * Math.Pow((x1 * y2 - x2 * y1 - x1 * y3 + x3 * y1 + x2 * y3 - x3 * y2), 2));
            ke[6, 2] = -(E * Math.Pow(t, 3) * ((2 * Math.Pow(x32, 2) - (2 * ga6 * x13 * x32) / my6) * (2 * Math.Pow(x13, 2) + nu * (2 * Math.Pow(y31, 2) - (2 * my5 * y23 * y31) / ga5) - (2 * my5 * x13 * x32) / ga5) + (2 * Math.Pow(y23, 2) - (2 * ga6 * y23 * y31) / my6) * (2 * Math.Pow(y31, 2) + nu * (2 * Math.Pow(x13, 2) - (2 * my5 * x13 * x32) / ga5) - (2 * my5 * y23 * y31) / ga5) - (nu / 2 - 1 / 2) * (4 * x13 * y31 - (my5 * (2 * x13 * y23 + 2 * x32 * y31)) / ga5) * (4 * x32 * y23 - (ga6 * (2 * x13 * y23 + 2 * x32 * y31)) / my6))) / (192 * Math.Pow(Area, 3) * (Math.Pow(nu, 2) - 1));
            ke[6, 6] = -(E * Math.Pow(t, 3) * ((2 * Math.Pow(x13, 2) - (2 * my5 * x13 * x32) / ga5) * (2 * Math.Pow(x13, 2) + nu * (2 * Math.Pow(y31, 2) - (2 * my5 * y23 * y31) / ga5) - (2 * my5 * x13 * x32) / ga5) - (nu / 2 - 1 / 2) * Math.Pow((4 * x13 * y31 - (my5 * (2 * x13 * y23 + 2 * x32 * y31)) / ga5), 2) + (2 * Math.Pow(y31, 2) - (2 * my5 * y23 * y31) / ga5) * (2 * Math.Pow(y31, 2) + nu * (2 * Math.Pow(x13, 2) - (2 * my5 * x13 * x32) / ga5) - (2 * my5 * y23 * y31) / ga5))) / (192 * Math.Pow(Area, 3) * (Math.Pow(nu, 2) - 1));
            ke[6, 7] = -(E * Math.Pow(t, 3) * ((4 * x13 * x32 * (ga5 * Math.Pow(x13, 2) + ga5 * nu * Math.Pow(y31, 2) - my5 * x13 * x32 - my5 * nu * y23 * y31)) / Math.Pow(ga5, 2) + (4 * y23 * y31 * (ga5 * Math.Pow(y31, 2) + ga5 * nu * Math.Pow(x13, 2) - my5 * y23 * y31 - my5 * nu * x13 * x32)) / Math.Pow(ga5, 2) + (2 * (2 * x13 * y23 + 2 * x32 * y31) * (nu / 2 - 1 / 2) * (my5 * x13 * y23 - 2 * ga5 * x13 * y31 + my5 * x32 * y31)) / Math.Pow(ga5, 2))) / (192 * Math.Pow(Area, 3) * (Math.Pow(nu, 2) - 1));
            ke[6, 10] = -(E * Math.Pow(t, 3) * (2 * x13 * x32 * (a5 / ga5 + a6 / my6) * (2 * Math.Pow(x13, 2) + nu * (2 * Math.Pow(y31, 2) - (2 * my5 * y23 * y31) / ga5) - (2 * my5 * x13 * x32) / ga5) - (2 * x13 * y23 + 2 * x32 * y31) * (nu / 2 - 1 / 2) * (4 * x13 * y31 - (my5 * (2 * x13 * y23 + 2 * x32 * y31)) / ga5) * (a5 / ga5 + a6 / my6) + 2 * y23 * y31 * (a5 / ga5 + a6 / my6) * (2 * Math.Pow(y31, 2) + nu * (2 * Math.Pow(x13, 2) - (2 * my5 * x13 * x32) / ga5) - (2 * my5 * y23 * y31) / ga5))) / (192 * Math.Pow(Area, 3) * (Math.Pow(nu, 2) - 1));
            ke[6, 11] = -(E * Math.Pow(t, 3) * ((2 * x13 * x32 * (2 * Math.Pow(x13, 2) + nu * (2 * Math.Pow(y31, 2) - (2 * my5 * y23 * y31) / ga5) - (2 * my5 * x13 * x32) / ga5)) / my6 + (2 * y23 * y31 * (2 * Math.Pow(y31, 2) + nu * (2 * Math.Pow(x13, 2) - (2 * my5 * x13 * x32) / ga5) - (2 * my5 * y23 * y31) / ga5)) / my6 - ((2 * x13 * y23 + 2 * x32 * y31) * (nu / 2 - 1 / 2) * (4 * x13 * y31 - (my5 * (2 * x13 * y23 + 2 * x32 * y31)) / ga5)) / my6)) / (192 * Math.Pow(Area, 3) * (Math.Pow(nu, 2) - 1));
            ke[7, 2] = (E * Math.Pow(t, 3) * ((4 * x32 * (x13 * x32 + nu * y23 * y31) * (ga6 * x13 - my6 * x32)) / (ga5 * my6) - (2 * (2 * x13 * y23 + 2 * x32 * y31) * (nu / 2 - 1 / 2) * (ga6 * x13 * y23 + ga6 * x32 * y31 - 2 * my6 * x32 * y23)) / (ga5 * my6) + (4 * y23 * (y23 * y31 + nu * x13 * x32) * (ga6 * y31 - my6 * y23)) / (ga5 * my6))) / (192 * Math.Pow(Area, 3) * (Math.Pow(nu, 2) - 1));
            ke[7, 6] = -(E * Math.Pow(t, 3) * ((2 * (2 * x13 * y23 + 2 * x32 * y31) * (nu / 2 - 1 / 2) * (my5 * x13 * y23 - 2 * ga5 * x13 * y31 + my5 * x32 * y31)) / Math.Pow(ga5, 2) + (4 * x13 * (x13 * x32 + nu * y23 * y31) * (ga5 * x13 - my5 * x32)) / Math.Pow(ga5, 2) + (4 * y31 * (y23 * y31 + nu * x13 * x32) * (ga5 * y31 - my5 * y23)) / Math.Pow(ga5, 2))) / (192 * Math.Pow(Area, 3) * (Math.Pow(nu, 2) - 1));
            ke[7, 7] = -(E * Math.Pow(t, 3) * ((2 * x13 * x32 * ((2 * x13 * x32) / ga5 + (2 * nu * y23 * y31) / ga5)) / ga5 - (Math.Pow((2 * x13 * y23 + 2 * x32 * y31), 2) * (nu / 2 - 1 / 2)) / Math.Pow(ga5, 2) + (2 * y23 * y31 * ((2 * y23 * y31) / ga5 + (2 * nu * x13 * x32) / ga5)) / ga5)) / (192 * Math.Pow(Area, 3) * (Math.Pow(nu, 2) - 1));
            ke[7, 10] = -(E * Math.Pow(t, 3) * (a6 * ga5 + a5 * my6) * (2 * Math.Pow(x13, 2) * Math.Pow(x32, 2) + Math.Pow(x13, 2) * Math.Pow(y23, 2) + Math.Pow(x32, 2) * Math.Pow(y31, 2) + 2 * Math.Pow(y23, 2) * Math.Pow(y31, 2) - nu * Math.Pow(x13, 2) * Math.Pow(y23, 2) - nu * Math.Pow(x32, 2) * Math.Pow(y31, 2) + 2 * x13 * x32 * y23 * y31 + 2 * nu * x13 * x32 * y23 * y31)) / (96 * Math.Pow(Area, 3) * Math.Pow(ga5, 2) * my6 * (Math.Pow(nu, 2) - 1));
            ke[7, 11] = -(E * Math.Pow(t, 3) * (2 * Math.Pow(x13, 2) * Math.Pow(x32, 2) + Math.Pow(x13, 2) * Math.Pow(y23, 2) + Math.Pow(x32, 2) * Math.Pow(y31, 2) + 2 * Math.Pow(y23, 2) * Math.Pow(y31, 2) - nu * Math.Pow(x13, 2) * Math.Pow(y23, 2) - nu * Math.Pow(x32, 2) * Math.Pow(y31, 2) + 2 * x13 * x32 * y23 * y31 + 2 * nu * x13 * x32 * y23 * y31)) / (96 * Math.Pow(Area, 3) * ga5 * my6 * (Math.Pow(nu, 2) - 1));
            ke[8, 0] = -(Area * E * t * (x1 * x2 - x1 * x3 + x2 * x3 + 2 * y1 * y2 - 2 * y1 * y3 + 2 * y2 * y3 + nu * Math.Pow(x2, 2) - Math.Pow(x2, 2) - 2 * Math.Pow(y2, 2) - nu * x1 * x2 + nu * x1 * x3 - nu * x2 * x3)) / (2 * (Math.Pow(nu, 2) - 1) * Math.Pow((x1 * y2 - x2 * y1 - x1 * y3 + x3 * y1 + x2 * y3 - x3 * y2), 2));
            ke[8, 1] = -(Area * E * t * (x1 * y3 - x1 * y2 + x2 * y2 - x2 * y3 + nu * x1 * y2 - 2 * nu * x2 * y1 - nu * x1 * y3 + nu * x2 * y2 + 2 * nu * x3 * y1 + nu * x2 * y3 - 2 * nu * x3 * y2)) / (2 * (Math.Pow(nu, 2) - 1) * Math.Pow((x1 * y2 - x2 * y1 - x1 * y3 + x3 * y1 + x2 * y3 - x3 * y2), 2));
            ke[8, 4] = -(Area * E * t * (x1 * x2 + x1 * x3 - x2 * x3 + 2 * y1 * y2 + 2 * y1 * y3 - 2 * y2 * y3 + nu * Math.Pow(x1, 2) - Math.Pow(x1, 2) - 2 * Math.Pow(y1, 2) - nu * x1 * x2 - nu * x1 * x3 + nu * x2 * x3)) / (2 * (Math.Pow(nu, 2) - 1) * Math.Pow((x1 * y2 - x2 * y1 - x1 * y3 + x3 * y1 + x2 * y3 - x3 * y2), 2));
            ke[8, 5] = -(Area * E * t * (x1 * y1 - x2 * y1 - x1 * y3 + x2 * y3 + nu * x1 * y1 - 2 * nu * x1 * y2 + nu * x2 * y1 + nu * x1 * y3 - 2 * nu * x3 * y1 - nu * x2 * y3 + 2 * nu * x3 * y2)) / (2 * (Math.Pow(nu, 2) - 1) * Math.Pow((x1 * y2 - x2 * y1 - x1 * y3 + x3 * y1 + x2 * y3 - x3 * y2), 2));
            ke[8, 8] = -(Area * E * t * (Math.Pow(x1, 2) - 4 * y1 * y2 - nu * Math.Pow(x1, 2) - nu * Math.Pow(x2, 2) - 2 * x1 * x2 + Math.Pow(x2, 2) + 2 * Math.Pow(y1, 2) + 2 * Math.Pow(y2, 2) + 2 * nu * x1 * x2)) / (2 * (Math.Pow(nu, 2) - 1) * Math.Pow((x1 * y2 - x2 * y1 - x1 * y3 + x3 * y1 + x2 * y3 - x3 * y2), 2));
            ke[8, 9] = (Area * E * t * (x1 - x2) * (y1 - y2) * (nu + 1)) / (2 * (Math.Pow(nu, 2) - 1) * Math.Pow((x1 * y2 - x2 * y1 - x1 * y3 + x3 * y1 + x2 * y3 - x3 * y2), 2));
            ke[9, 0] = -(Area * E * t * (x2 * y2 - x2 * y1 + x3 * y1 - x3 * y2 - 2 * nu * x1 * y2 + nu * x2 * y1 + 2 * nu * x1 * y3 + nu * x2 * y2 - nu * x3 * y1 - 2 * nu * x2 * y3 + nu * x3 * y2)) / (2 * (Math.Pow(nu, 2) - 1) * Math.Pow((x1 * y2 - x2 * y1 - x1 * y3 + x3 * y1 + x2 * y3 - x3 * y2), 2));
            ke[9, 1] = -(Area * E * t * (2 * x1 * x2 - 2 * x1 * x3 + 2 * x2 * x3 + y1 * y2 - y1 * y3 + y2 * y3 + nu * Math.Pow(y2, 2) - 2 * Math.Pow(x2, 2) - Math.Pow(y2, 2) - nu * y1 * y2 + nu * y1 * y3 - nu * y2 * y3)) / (2 * (Math.Pow(nu, 2) - 1) * Math.Pow((x1 * y2 - x2 * y1 - x1 * y3 + x3 * y1 + x2 * y3 - x3 * y2), 2));
            ke[9, 4] = -(Area * E * t * (x1 * y1 - x1 * y2 - x3 * y1 + x3 * y2 + nu * x1 * y1 + nu * x1 * y2 - 2 * nu * x2 * y1 - 2 * nu * x1 * y3 + nu * x3 * y1 + 2 * nu * x2 * y3 - nu * x3 * y2)) / (2 * (Math.Pow(nu, 2) - 1) * Math.Pow((x1 * y2 - x2 * y1 - x1 * y3 + x3 * y1 + x2 * y3 - x3 * y2), 2));
            ke[9, 5] = -(Area * E * t * (2 * x1 * x2 + 2 * x1 * x3 - 2 * x2 * x3 + y1 * y2 + y1 * y3 - y2 * y3 + nu * Math.Pow(y1, 2) - 2 * Math.Pow(x1, 2) - Math.Pow(y1, 2) - nu * y1 * y2 - nu * y1 * y3 + nu * y2 * y3)) / (2 * (Math.Pow(nu, 2) - 1) * Math.Pow((x1 * y2 - x2 * y1 - x1 * y3 + x3 * y1 + x2 * y3 - x3 * y2), 2));
            ke[9, 8] = (Area * E * t * (x1 - x2) * (y1 - y2) * (nu + 1)) / (2 * (Math.Pow(nu, 2) - 1) * Math.Pow((x1 * y2 - x2 * y1 - x1 * y3 + x3 * y1 + x2 * y3 - x3 * y2), 2));
            ke[9, 9] = -(Area * E * t * (2 * Math.Pow(x1, 2) - 2 * y1 * y2 - nu * Math.Pow(y1, 2) - nu * Math.Pow(y2, 2) - 4 * x1 * x2 + 2 * Math.Pow(x2, 2) + Math.Pow(y1, 2) + Math.Pow(y2, 2) + 2 * nu * y1 * y2)) / (2 * (Math.Pow(nu, 2) - 1) * Math.Pow((x1 * y2 - x2 * y1 - x1 * y3 + x3 * y1 + x2 * y3 - x3 * y2), 2));
            ke[10, 2] = -(E * Math.Pow(t, 3) * ((2 * Math.Pow(x32, 2) - (2 * ga6 * x13 * x32) / my6) * (2 * x13 * x32 * (a5 / ga5 + a6 / my6) + 2 * nu * y23 * y31 * (a5 / ga5 + a6 / my6)) + (2 * Math.Pow(y23, 2) - (2 * ga6 * y23 * y31) / my6) * (2 * y23 * y31 * (a5 / ga5 + a6 / my6) + 2 * nu * x13 * x32 * (a5 / ga5 + a6 / my6)) - (2 * x13 * y23 + 2 * x32 * y31) * (nu / 2 - 1 / 2) * (4 * x32 * y23 - (ga6 * (2 * x13 * y23 + 2 * x32 * y31)) / my6) * (a5 / ga5 + a6 / my6))) / (192 * Math.Pow(Area, 3) * (Math.Pow(nu, 2) - 1));
            ke[10, 6] = -(E * Math.Pow(t, 3) * ((2 * Math.Pow(x13, 2) - (2 * my5 * x13 * x32) / ga5) * (2 * x13 * x32 * (a5 / ga5 + a6 / my6) + 2 * nu * y23 * y31 * (a5 / ga5 + a6 / my6)) + (2 * Math.Pow(y31, 2) - (2 * my5 * y23 * y31) / ga5) * (2 * y23 * y31 * (a5 / ga5 + a6 / my6) + 2 * nu * x13 * x32 * (a5 / ga5 + a6 / my6)) - (2 * x13 * y23 + 2 * x32 * y31) * (nu / 2 - 1 / 2) * (4 * x13 * y31 - (my5 * (2 * x13 * y23 + 2 * x32 * y31)) / ga5) * (a5 / ga5 + a6 / my6))) / (192 * Math.Pow(Area, 3) * (Math.Pow(nu, 2) - 1));
            ke[10, 7] = -(E * Math.Pow(t, 3) * (a6 * ga5 + a5 * my6) * (2 * Math.Pow(x13, 2) * Math.Pow(x32, 2) + Math.Pow(x13, 2) * Math.Pow(y23, 2) + Math.Pow(x32, 2) * Math.Pow(y31, 2) + 2 * Math.Pow(y23, 2) * Math.Pow(y31, 2) - nu * Math.Pow(x13, 2) * Math.Pow(y23, 2) - nu * Math.Pow(x32, 2) * Math.Pow(y31, 2) + 2 * x13 * x32 * y23 * y31 + 2 * nu * x13 * x32 * y23 * y31)) / (96 * Math.Pow(Area, 3) * Math.Pow(ga5, 2) * my6 * (Math.Pow(nu, 2) - 1));
            ke[10, 10] = -(E * Math.Pow(t, 3) * Math.Pow((a6 * ga5 + a5 * my6), 2) * (2 * Math.Pow(x13, 2) * Math.Pow(x32, 2) + Math.Pow(x13, 2) * Math.Pow(y23, 2) + Math.Pow(x32, 2) * Math.Pow(y31, 2) + 2 * Math.Pow(y23, 2) * Math.Pow(y31, 2) - nu * Math.Pow(x13, 2) * Math.Pow(y23, 2) - nu * Math.Pow(x32, 2) * Math.Pow(y31, 2) + 2 * x13 * x32 * y23 * y31 + 2 * nu * x13 * x32 * y23 * y31)) / (96 * Math.Pow(Area, 3) * Math.Pow(ga5, 2) * Math.Pow(my6, 2) * (Math.Pow(nu, 2) - 1));
            ke[10, 11] = -(E * Math.Pow(t, 3) * (a6 * ga5 + a5 * my6) * (2 * Math.Pow(x13, 2) * Math.Pow(x32, 2) + Math.Pow(x13, 2) * Math.Pow(y23, 2) + Math.Pow(x32, 2) * Math.Pow(y31, 2) + 2 * Math.Pow(y23, 2) * Math.Pow(y31, 2) - nu * Math.Pow(x13, 2) * Math.Pow(y23, 2) - nu * Math.Pow(x32, 2) * Math.Pow(y31, 2) + 2 * x13 * x32 * y23 * y31 + 2 * nu * x13 * x32 * y23 * y31)) / (96 * Math.Pow(Area, 3) * ga5 * Math.Pow(my6, 2) * (Math.Pow(nu, 2) - 1));
            ke[11, 2] = (E * Math.Pow(t, 3) * ((4 * x32 * (x13 * x32 + nu * y23 * y31) * (ga6 * x13 - my6 * x32)) / Math.Pow(my6, 2) - (2 * (2 * x13 * y23 + 2 * x32 * y31) * (nu / 2 - 1 / 2) * (ga6 * x13 * y23 + ga6 * x32 * y31 - 2 * my6 * x32 * y23)) / Math.Pow(my6, 2) + (4 * y23 * (y23 * y31 + nu * x13 * x32) * (ga6 * y31 - my6 * y23)) / Math.Pow(my6, 2))) / (192 * Math.Pow(Area, 3) * (Math.Pow(nu, 2) - 1));
            ke[11, 6] = -(E * Math.Pow(t, 3) * ((2 * (2 * x13 * y23 + 2 * x32 * y31) * (nu / 2 - 1 / 2) * (my5 * x13 * y23 - 2 * ga5 * x13 * y31 + my5 * x32 * y31)) / (ga5 * my6) + (4 * x13 * (x13 * x32 + nu * y23 * y31) * (ga5 * x13 - my5 * x32)) / (ga5 * my6) + (4 * y31 * (y23 * y31 + nu * x13 * x32) * (ga5 * y31 - my5 * y23)) / (ga5 * my6))) / (192 * Math.Pow(Area, 3) * (Math.Pow(nu, 2) - 1));
            ke[11, 7] = -(E * Math.Pow(t, 3) * (2 * Math.Pow(x13, 2) * Math.Pow(x32, 2) + Math.Pow(x13, 2) * Math.Pow(y23, 2) + Math.Pow(x32, 2) * Math.Pow(y31, 2) + 2 * Math.Pow(y23, 2) * Math.Pow(y31, 2) - nu * Math.Pow(x13, 2) * Math.Pow(y23, 2) - nu * Math.Pow(x32, 2) * Math.Pow(y31, 2) + 2 * x13 * x32 * y23 * y31 + 2 * nu * x13 * x32 * y23 * y31)) / (96 * Math.Pow(Area, 3) * ga5 * my6 * (Math.Pow(nu, 2) - 1));
            ke[11, 10] = -(E * Math.Pow(t, 3) * (a6 * ga5 + a5 * my6) * (2 * Math.Pow(x13, 2) * Math.Pow(x32, 2) + Math.Pow(x13, 2) * Math.Pow(y23, 2) + Math.Pow(x32, 2) * Math.Pow(y31, 2) + 2 * Math.Pow(y23, 2) * Math.Pow(y31, 2) - nu * Math.Pow(x13, 2) * Math.Pow(y23, 2) - nu * Math.Pow(x32, 2) * Math.Pow(y31, 2) + 2 * x13 * x32 * y23 * y31 + 2 * nu * x13 * x32 * y23 * y31)) / (96 * Math.Pow(Area, 3) * ga5 * Math.Pow(my6, 2) * (Math.Pow(nu, 2) - 1));
            ke[11, 11] = -(E * Math.Pow(t, 3) * ((2 * x13 * x32 * ((2 * x13 * x32) / my6 + (2 * nu * y23 * y31) / my6)) / my6 - (Math.Pow((2 * x13 * y23 + 2 * x32 * y31), 2) * (nu / 2 - 1 / 2)) / Math.Pow(my6, 2) + (2 * y23 * y31 * ((2 * y23 * y31) / my6 + (2 * nu * x13 * x32) / my6)) / my6)) / (192 * Math.Pow(Area, 3) * (Math.Pow(nu, 2) - 1));

            Matrix<double> Ke = ke.Multiply(T);
            Ke = T_T.Multiply(Ke);
             
            return Ke;
        }


        private List<double> CreateLoadList(List<string> loadtxt, List<string> momenttxt, List<Point3d> vertices)
        {
            List<double> loads = new List<double>(new double[vertices.Count * 6]);
            List<double> inputLoads = new List<double>();
            List<Point3d> coordlist = new List<Point3d>();

            for (int i = 0; i < loadtxt.Count; i++)
            {
                string coordstr = (loadtxt[i].Split(':')[0]);
                string loadstr = (loadtxt[i].Split(':')[1]);

                string[] coordstr1 = (coordstr.Split(','));
                string[] loadstr1 = (loadstr.Split(','));

                inputLoads.Add(Math.Round(double.Parse(loadstr1[0]), 2));
                inputLoads.Add(Math.Round(double.Parse(loadstr1[1]), 2));
                inputLoads.Add(Math.Round(double.Parse(loadstr1[2]), 2));

                coordlist.Add(new Point3d(Math.Round(double.Parse(coordstr1[0]), 2), Math.Round(double.Parse(coordstr1[1]), 2), Math.Round(double.Parse(coordstr1[2]), 2)));
            }

            foreach (Point3d point in coordlist)
            {
                int i = vertices.IndexOf(point);
                int j = coordlist.IndexOf(point);
                loads[i * 6 + 0] = inputLoads[j * 3 + 0];
                loads[i * 6 + 1] = inputLoads[j * 3 + 1];
                loads[i * 6 + 2] = inputLoads[j * 3 + 2];
            }
            inputLoads.Clear();
            coordlist.Clear();
            for (int i = 0; i < momenttxt.Count; i++) if (momenttxt[0] != "")
                {
                    string coordstr = (momenttxt[i].Split(':')[0]);
                    string loadstr = (momenttxt[i].Split(':')[1]);

                    string[] coordstr1 = (coordstr.Split(','));
                    string[] loadstr1 = (loadstr.Split(','));

                    inputLoads.Add(Math.Round(double.Parse(loadstr1[0]), 2));
                    inputLoads.Add(Math.Round(double.Parse(loadstr1[1]), 2));
                    inputLoads.Add(Math.Round(double.Parse(loadstr1[2]), 2));


                    coordlist.Add(new Point3d(Math.Round(double.Parse(coordstr1[0]), 2), Math.Round(double.Parse(coordstr1[1]), 2), Math.Round(double.Parse(coordstr1[2]), 2)));
                }

            foreach (Point3d point in coordlist)
            {
                int i = vertices.IndexOf(point);
                int j = coordlist.IndexOf(point);
                loads[i * 6 + 3] = inputLoads[j * 3 + 0];
                loads[i * 6 + 4] = inputLoads[j * 3 + 1];
                loads[i * 6 + 5] = inputLoads[j * 3 + 2];
            }
            return loads;
        }

        private List<int> CreateBDCList(List<string> bdctxt, List<Point3d> vertices)
        {
            List<int> bdc_value = new List<int>(new int[vertices.Count * 6]);
            List<int> bdcs = new List<int>();
            List<Point3d> bdc_points = new List<Point3d>(); //Coordinates relating til bdc_value in for (eg. x y z)

            //Parse string input
            for (int i = 0; i < bdctxt.Count; i++)
            {
                string coordstr = (bdctxt[i].Split(':')[0]);
                string bdcstr = (bdctxt[i].Split(':')[1]);

                string[] coordstr1 = (coordstr.Split(','));
                string[] bdcstr1 = (bdcstr.Split(','));

                bdc_points.Add(new Point3d(Math.Round(double.Parse(coordstr1[0]), 2), Math.Round(double.Parse(coordstr1[1]), 2), Math.Round(double.Parse(coordstr1[2]), 2)));

                bdcs.Add(int.Parse(bdcstr1[0]));
                bdcs.Add(int.Parse(bdcstr1[1]));
                bdcs.Add(int.Parse(bdcstr1[2]));
                bdcs.Add(int.Parse(bdcstr1[3]));
                bdcs.Add(int.Parse(bdcstr1[4]));
                bdcs.Add(int.Parse(bdcstr1[5]));
            }


            //Format to correct entries in bdc_value
            for (int i = 0; i < vertices.Count; i++)
            {
                Point3d tempP = vertices[i];

                if (bdc_points.Contains(tempP))
                {
                    bdc_value[i * 6 + 0] = bdcs[bdc_points.IndexOf(tempP) * 6 + 0];
                    bdc_value[i * 6 + 1] = bdcs[bdc_points.IndexOf(tempP) * 6 + 1];
                    bdc_value[i * 6 + 2] = bdcs[bdc_points.IndexOf(tempP) * 6 + 2];
                    bdc_value[i * 6 + 3] = bdcs[bdc_points.IndexOf(tempP) * 6 + 3];
                    bdc_value[i * 6 + 4] = bdcs[bdc_points.IndexOf(tempP) * 6 + 4];
                    bdc_value[i * 6 + 5] = bdcs[bdc_points.IndexOf(tempP) * 6 + 5];
                }
                else
                {
                    bdc_value[i * 6 + 0] = 1;
                    bdc_value[i * 6 + 1] = 1;
                    bdc_value[i * 6 + 2] = 1;
                    bdc_value[i * 6 + 3] = 1;
                    bdc_value[i * 6 + 4] = 1;
                    bdc_value[i * 6 + 5] = 1;
                }
            }
            return bdc_value;
        }

        private void SetMaterial(string mattxt, out double E, out double A, out double Iy, out double Iz, out double J, out double G, out double nu, out double t)
        {
            string[] matProp = (mattxt.Split(','));

            E = (Math.Round(double.Parse(matProp[0]), 2));
            A = (Math.Round(double.Parse(matProp[1]), 2));
            Iy = (Math.Round(double.Parse(matProp[2]), 2));
            Iz = (Math.Round(double.Parse(matProp[3]), 2));
            G = (Math.Round(double.Parse(matProp[4]), 2));
            nu = (Math.Round(double.Parse(matProp[5]), 3));
            t = (Math.Round(double.Parse(matProp[6]), 2));
            J = Iy + Iz;
        }

        protected override System.Drawing.Bitmap Icon
        {
            get
            {

                return null;
            }
        }

        public override Guid ComponentGuid
        {
            get { return new Guid("3a61d696-911f-46cd-a687-ef48a48575b0"); }
        }
    }
}
