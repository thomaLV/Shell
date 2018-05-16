using System;
using System.Collections.Generic;

using Grasshopper.Kernel;
using Rhino.Geometry;
using System.Drawing;
using Grasshopper.GUI.Canvas;
using System.Windows.Forms;
using Grasshopper.GUI;

using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;
using System.Diagnostics;

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

        static bool startCalc = false;

        public static void setStart(string s, bool i)
        {
            if (s == "Run")
            {
                startCalc = i;
            }
        }

        public override void CreateAttributes()
        {
            m_attributes = new Attributes_Custom(this);
        }

        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddMeshParameter("IsoMesh", "IM", "The shell 6-node element isoparametric mesh, made by IsoMesher component", GH_ParamAccess.item);
            pManager.AddTextParameter("Boundary Conditions", "BDC", "Boundary Conditions in form x,y,z,vx,vy,vz,rx,ry,rz", GH_ParamAccess.list);
            pManager.AddTextParameter("Material properties", "Mat", "Material Properties", GH_ParamAccess.item, "210000,3600,4920000,4920000,79300,0.3,10");
            pManager.AddTextParameter("PointLoads", "PL", "Load given as Vector [N]", GH_ParamAccess.list);
        }

        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddNumberParameter("Deformations", "Def", "Deformations", GH_ParamAccess.list);
            pManager.AddTextParameter("Reactions", "R", "Reaction Forces", GH_ParamAccess.item);
            pManager.AddTextParameter("Element stresses", "Strs", "The Stress in each element", GH_ParamAccess.list);
            pManager.AddTextParameter("Element strains", "Strn", "The Strain in each element", GH_ParamAccess.list);
            pManager.AddLineParameter("Edges", "edges", "", GH_ParamAccess.list);
            pManager.AddNumberParameter("naked edge index", "", "", GH_ParamAccess.list);
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

            if (!DA.GetData(0, ref mesh)) return;           //sets inputted mesh into variable
            if (!DA.GetDataList(1, bdctxt)) return;         //sets boundary conditions as string
            if (!DA.GetData(2, ref mattxt)) return;         //sets material properties as string
            if (!DA.GetDataList(3, loadtxt)) return;        //sets load as string

            foreach (var face in mesh.Faces)
            {
                faces.Add(face);
            }

            foreach (var vertice in mesh.Vertices)
            {
                Point3d temp_vertice = new Point3d();
                temp_vertice.X = Math.Round(vertice.X, 4);
                temp_vertice.Y = Math.Round(vertice.Y, 4);
                temp_vertice.Z = Math.Round(vertice.Z, 4);
                vertices.Add(temp_vertice);
            }

            // Number of edges from Euler's formula
            int NoOfEdges = vertices.Count + faces.Count - 1;
            List<Line> edges = new List<Line>(NoOfEdges);
            Vector<double> nakedEdge = Vector<double>.Build.Dense(NoOfEdges,1);
            foreach (var face in faces)
            {
                Point3d vA = vertices[face.A];
                Point3d vB = vertices[face.B];
                Point3d vC = vertices[face.C];
                Line lineAB = new Line(vA, vB);
                Line lineBA = new Line(vB, vA);
                Line lineCB = new Line(vC, vB);
                Line lineBC = new Line(vB, vC);
                Line lineAC = new Line(vA, vC);
                Line lineCA = new Line(vC, vA);

                if (!edges.Contains(lineAB) && !edges.Contains(lineBA))
                {
                    edges.Add(lineAB);
                }
                else
                {
                    int i = edges.IndexOf(lineAB);
                    if (i == -1)
                    {
                        i = edges.IndexOf(lineBA);
                    }
                    nakedEdge[i] = 0;
                }
                if (!edges.Contains(lineCB) && !edges.Contains(lineBC))
                {
                    edges.Add(lineBC);
                }
                else
                {
                    int i = edges.IndexOf(lineBC);
                    if (i == -1)
                    {
                        i = edges.IndexOf(lineCB);
                    }
                    nakedEdge[i] = 0;
                }
                if (!edges.Contains(lineAC) && !edges.Contains(lineCA))
                {
                    edges.Add(lineAC);
                }
                else
                {
                    int i = edges.IndexOf(lineAC);
                    if (i == -1)
                    {
                        i = edges.IndexOf(lineCA);
                    }
                    nakedEdge[i] = 0;
                }
            }

            List<Point3d> uniqueNodes;
            GetUniqueNodes(vertices, out uniqueNodes);
            int gdofs = uniqueNodes.Count * 3 + edges.Count;

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

            Vector<double> def_tot;
            Vector<double> reactions;
            Vector<double> internalStresses;
            Vector<double> internalStrains;

            #region Prepares boundary conditions and loads for calculation

            //Interpret the BDC inputs (text) and create list of boundary condition (1/0 = free/clamped) for each dof.
            Vector<double> bdc_value = CreateBDCList(bdctxt, uniqueNodes, faces, vertices, edges);

            Vector<double> nakededge = Vector<double>.Build.Dense(gdofs, 0);
            for (int i = uniqueNodes.Count*3; i < gdofs; i++)
            {
                if (bdc_value[i] == 1)
                {
                    nakededge[i] = (nakedEdge[i - uniqueNodes.Count * 3]);
                }
            }
            List<double> test1 = new List<double>(nakededge.ToArray());

            //Interpreting input load (text) and creating load list (double)
            List<double> load = CreateLoadList(loadtxt, momenttxt, uniqueNodes, faces, vertices, edges);
            #endregion

            Matrix<double> K_red;
            Vector<double> load_red;

            String time = "TrySolve Start:" + Environment.NewLine;
            long timer = 0;
            Stopwatch watch = new Stopwatch();


            #region Create global and reduced stiffness matrix

            //Create global stiffness matrix

            watch.Start();
            Matrix<double> K_tot = GlobalStiffnessMatrix(faces, vertices, edges, uniqueNodes, E, A, Iy, Iz, J, G, nu, t);
            watch.Stop();
            timer = watch.ElapsedMilliseconds - timer;
            time += "Global stiffness matrix assembly: " + timer.ToString() + Environment.NewLine;

            //Create reduced K-matrix and reduced load list (removed clamped dofs)

            watch.Start();
            CreateReducedGlobalStiffnessMatrix(bdc_value, K_tot, load, uniqueNodes, nakededge, out K_red, out load_red);
            watch.Stop();
            timer = watch.ElapsedMilliseconds - timer;
            time += "Reduce global stiffness matrix: " + timer.ToString() + Environment.NewLine;
            
            #endregion

            if (startCalc)
            {

                #region Calculate deformations, reaction forces and internal strains and stresses

                bool test = K_red.IsSymmetric();
                //Calculate deformations
                Vector<double> def_reduced = Vector<double>.Build.Dense(K_red.ColumnCount);
                    watch.Start();
                    //def_reduced = K_red.Cholesky().Solve(load_red);
                def_reduced = K_red.Solve(load_red);
                    watch.Stop();
                    timer = watch.ElapsedMilliseconds - timer;
                    time += "Cholesky solve: " + timer.ToString() + Environment.NewLine;

                //Add the clamped dofs (= 0) to the deformations list
                def_tot = RestoreTotalDeformationVector(def_reduced, bdc_value, nakededge);

                //Calculate the reaction forces from the deformations
                reactions = K_tot.Multiply(def_tot);

                //Calculate the internal strains and stresses in each member
                // m = -h^3/12 * C * Bk * v
                // 
                // strain = B * v
                // stress = C * strain

                //CalculateInternalStrainsAndStresses(def_tot, vertices, E, out internalStresses, out internalStrains);
                #endregion
            }
            else
            {
                def_tot = Vector<double>.Build.Dense(bdc_value.Count * 6);
                //reactions = def_tot;

                //internalStresses = Vector<double>.Build.Dense(bdc_value.Count * 6);
                //internalStrains = internalStresses;
            }


            DA.SetDataList(0, def_tot);
            DA.SetData(1, time.ToString());
            DA.SetData(2, K_red.ToString(32,32));
            DA.SetData(3, K_tot.ToString(43,43));
            DA.SetDataList(4, edges);           
            DA.SetDataList(5, test1);
        }

        //private void CalculateInternalStrainsAndStresses(Vector<double> def, List<Point3d> vertices, double E, out Vector<double> internalStresses, out Vector<double> internalStrains)
        //{
        //    //preallocating lists
        //    internalStresses = new List<double>(geometry.Count);
        //    internalStrains = new List<double>(geometry.Count);

        //    foreach (Line line in geometry)
        //    {
        //        int index1 = points.IndexOf(new Point3d(Math.Round(line.From.X, 5), Math.Round(line.From.Y, 5), Math.Round(line.From.Z, 5)));
        //        int index2 = points.IndexOf(new Point3d(Math.Round(line.To.X, 5), Math.Round(line.To.Y, 5), Math.Round(line.To.Z, 5)));

        //        //fetching deformation of point
        //        double x1 = def[index1 * 3 + 0];
        //        double y1 = def[index1 * 3 + 1];
        //        double z1 = def[index1 * 3 + 2];
        //        double x2 = def[index2 * 3 + 0];
        //        double y2 = def[index2 * 3 + 1];
        //        double z2 = def[index2 * 3 + 2];

        //        //new node coordinates for deformed nodes
        //        double nx1 = points[index1].X + x1;
        //        double ny1 = points[index1].X + y1;
        //        double nz1 = points[index1].Z + z1;
        //        double nx2 = points[index2].X + x2;
        //        double ny2 = points[index2].X + y2;
        //        double nz2 = points[index2].Z + z2;

        //        //calculating dL = length of deformed line - original length of line
        //        double dL = Math.Sqrt(Math.Pow((nx2 - nx1), 2) + Math.Pow((ny2 - ny1), 2) + Math.Pow((nz2 - nz1), 2)) - line.Length;

        //        //calculating strain and stress
        //        internalStrains.Add(dL / line.Length);
        //        internalStresses.Add(internalStrains[internalStrains.Count - 1] * E);
        //    }
        //}

        private Vector<double> RestoreTotalDeformationVector(Vector<double> deformations_red, Vector<double> bdc_value, Vector<double> nakededges)
        {
            Vector<double> def = Vector<double>.Build.Dense(bdc_value.Count);
            for (int i = 0, j = 0; i < bdc_value.Count; i++)
            {
                if (bdc_value[i] == 1 && nakededges[i] == 0)
                {
                    def[i] = deformations_red[j];
                    j++;
                }
            }
            return def;
        }

        private void CreateReducedGlobalStiffnessMatrix(Vector<double> bdc_value, Matrix<double> K, List<double> load, List<Point3d> uniqueNodes, Vector<double> nakededges, out Matrix<double> K_red, out Vector<double> load_red)
        {
            int oldRC = load.Count;
            int newRC = Convert.ToInt16(bdc_value.Sum()-nakededges.Sum());
            K_red = Matrix<double>.Build.Dense(newRC, newRC, 0);
            load_red = Vector<double>.Build.Dense(newRC, 0);
            for (int i = 0, ii = 0; i < oldRC; i++)
            {
                //is bdc_value in row i free?
                if (bdc_value[i] == 1 && nakededges[i] == 0)
                {                    
                    for (int j = 0, jj = 0; j < oldRC; j++)
                    {
                        //is bdc_value in col j free?
                        if (bdc_value[j] == 1 && nakededges[j] == 0)
                        {                                
                            //if yes, then add to new K
                            K_red[i - ii, j - jj] = Math.Round(K[i, j], 4);                               
                        }
                        //else if (bdc_value[j] == 0 && nakededges[j] == 1)
                        //{
                        //    //if not free, remember to skip 1 column when adding next time
                        //    jj++;
                        //    jj++;
                        //}
                        else
                        {
                            jj++;
                        }
                    }
                    //add to reduced load list
                    load_red[i - ii] = load[i];
                }
                //else if (bdc_value[i] == 0 && nakededges[i] == 1)
                //{
                //    //if not free, remember to skip 1 row when adding next time
                //    ii++;
                //    ii++;
                //}
                else
                {
                    ii++;
                }
            }
            //for (int i = 0, j=0; i < size; i++)
            //{
            //    //remove clamped dofs
            //    if (bdc_value[i] == 0)
            //    {
            //        K_red = K_red.RemoveRow(i - j);
            //        K_red = K_red.RemoveColumn(i - j);
            //        load_redu.RemoveAt(i - j);
            //        j++;
            //    }
            //}
            //load_red = Vector<double>.Build.DenseOfEnumerable(load_redu);
        }

        private void GetUniqueNodes(List<Point3d> vertices, out List<Point3d> uniqueNodes)
        {
            uniqueNodes = new List<Point3d>();
            for (int i = 0; i < vertices.Count; i++)
            {
                Point3d tempNode = new Point3d(Math.Round(vertices[i].X, 4), Math.Round(vertices[i].Y, 4), Math.Round(vertices[i].Z, 4));
                if (!uniqueNodes.Contains(tempNode))
                {
                    uniqueNodes.Add(tempNode);
                }
            }
        }

        private Matrix<double> GlobalStiffnessMatrix(List<MeshFace> faces, List<Point3d> vertices, List<Line> edges, List<Point3d> uniqueNodes, double E, double A, double Iy, double Iz, double J, double G, double nu, double t)
        {
            int gdofs = uniqueNodes.Count * 3 + edges.Count;

            // take into account that there are more sides than veritces

            var KG = Matrix<double>.Build.Dense(gdofs, gdofs);

            foreach (var face in faces)
            {
                int indexA = uniqueNodes.IndexOf(vertices[face.A]);
                int indexB = uniqueNodes.IndexOf(vertices[face.B]); 
                int indexC = uniqueNodes.IndexOf(vertices[face.C]);

                Point3d verticeA = vertices[indexA];
                Point3d verticeB = vertices[indexB];
                Point3d verticeC = vertices[indexC];

                int edgeIndex1 = edges.IndexOf(new Line(verticeA, verticeB));
                if (edgeIndex1 == -1) { edgeIndex1 = edges.IndexOf(new Line(verticeB, verticeA)); }
                int edgeIndex2 = edges.IndexOf(new Line(verticeB, verticeC));
                if (edgeIndex2 == -1) { edgeIndex2 = edges.IndexOf(new Line(verticeC, verticeB)); }
                int edgeIndex3 = edges.IndexOf(new Line(verticeC, verticeA));
                if (edgeIndex3 == -1) { edgeIndex3 = edges.IndexOf(new Line(verticeA, verticeC)); }

                int[] eindx = new int[] { edgeIndex1, edgeIndex2, edgeIndex3 };
                int[] vindx = new int[] { indexA, indexB, indexC };

                double x1 = verticeA.X;
                double x2 = verticeB.X;
                double x3 = verticeC.X;

                double y1 = verticeA.Y;
                double y2 = verticeB.Y;
                double y3 = verticeC.Y;

                double z1 = verticeA.Z;
                double z2 = verticeB.Z;
                double z3 = verticeC.Z;

                double[] xList = new double[3] { x1, x2, x3 };
                double[] yList = new double[3] { y1, y2, y3 };
                double[] zList = new double[3] { z1, z2, z3 };

                Matrix<double> Ke = ElementStiffnessMatrix(xList, yList, zList, E, nu, t);

                int nodeDofs = uniqueNodes.Count * 3;
                for (int row = 0; row < 3; row++)
                {
                    for (int col = 0; col < 3; col++)
                    {
                        //top left 3x3 of K-element matrix
                        KG[indexA * 3 + row, indexA * 3 + col] += Ke[row, col];
                        //top middle 3x3 of k-element matrix
                        KG[indexA * 3 + row, indexB * 3 + col] += Ke[row, col + 4];
                        //top right 3x3 of k-element matrix 
                        KG[indexA * 3 + row, indexC * 3 + col] += Ke[row, col + 4 * 2];

                        //middle left 3x3 of k-element matrix
                        KG[indexB * 3 + row, indexA * 3 + col] += Ke[row + 4, col];
                        //middle middle 3x3 of k-element matrix
                        KG[indexB * 3 + row, indexB * 3 + col] += Ke[row + 4, col + 4];
                        //middle right 3x3 of k-element matrix
                        KG[indexB * 3 + row, indexC * 3 + col] += Ke[row + 4, col + 4 * 2];

                        //bottom left 3x3 of k-element matrix
                        KG[indexC * 3 + row, indexA * 3 + col] += Ke[row + 4 * 2, col];
                        //bottom middle 3x3 of k-element matrix
                        KG[indexC * 3 + row, indexB * 3 + col] += Ke[row + 4 * 2, col + 4];
                        //bottom right 3x3 of k-element matrix
                        KG[indexC * 3 + row, indexC * 3 + col] += Ke[row + 4 * 2, col + 4 * 2];

                        // insert rotations for edges in correct place
                        //Rotation to rotation relation
                        //double insrt = Ke[row * 4 + 3, col * 4 + 3];
                        //if (insrt == 0)
                        //{
                        //    throw new System.ArgumentException("digonal cannot be zero");
                        //}
                        KG[nodeDofs + eindx[row], nodeDofs + eindx[col]] += Ke[row * 4 + 3, col * 4 + 3];
                        //Rotation to z relation lower left
                        KG[nodeDofs + eindx[row], vindx[col] * 3 + 2] += Ke[row * 4 + 3, col * 4 + 2];
                        //Rotation to z relation lower left
                        KG[vindx[row] * 3 + 2, nodeDofs + eindx[col]] += Ke[row * 4 + 2, col * 4 + 3];
                    }
                }



                    //int ldofs = 4;

                    //Inputting values to correct entries in Global Stiffness Matrix
                    //for (int row = 0; row < ldofs; row++)
                    //{
                    //    for (int col = 0; col < ldofs; col++)
                    //    {
                    //        //top left 4x4 of K-element matrix
                    //        KG[indexA * ldofs + row, indexA * ldofs + col] += Ke[row, col];
                    //        //top middle 4x4 of k-element matrix
                    //        KG[indexA * ldofs + row, indexB * ldofs + col] += Ke[row, col + ldofs];
                    //        //top right 4x4 of k-element matrix  
                    //        KG[indexA * ldofs + row, indexC * ldofs + col] += Ke[row, col + ldofs * 2];

                    //        //middle left 4x4 of k-element matrix
                    //        KG[indexB * ldofs + row, indexA * ldofs + col] += Ke[row + ldofs, col];
                    //        //middle middle 4x4 of k-element matrix
                    //        KG[indexB * ldofs + row, indexB * ldofs + col] += Ke[row + ldofs, col + ldofs];
                    //        //middle right 4x4 of k-element matrix
                    //        KG[indexB * ldofs + row, indexC * ldofs + col] += Ke[row + ldofs, col + ldofs * 2];

                    //        //bottom left 4x4 of k-element matrix
                    //        KG[indexC * ldofs + row, indexA * ldofs + col] += Ke[row + ldofs * 2, col];
                    //        //bottom middle 4x4 of k-element matrix
                    //        KG[indexC * ldofs + row, indexB * ldofs + col] += Ke[row + ldofs * 2, col + ldofs];
                    //        //bottom right 4x4 of k-element matrix
                    //        KG[indexC * ldofs + row, indexC * ldofs + col] += Ke[row + ldofs * 2, col + ldofs * 2];
                    //    }
                    //}
                }
            return KG;
        }

        private Matrix<double> ElementStiffnessMatrix(double[] xList, double[] yList, double[] zList, double E, double nu, double t)
        {

            #region Get global coordinates and transform into local cartesian system

            // fetching global coordinates
            double x1 = xList[0];
            double x2 = xList[1];
            double x3 = xList[2];

            double y1 = yList[0];
            double y2 = yList[1];
            double y3 = yList[2];

            double z1 = zList[0];
            double z2 = zList[1];
            double z3 = zList[2];

            // determine angles for tranformation matrix
            double Lx = Math.Sqrt((Math.Pow((x1 - x2), 2) + Math.Pow((y1 - y2), 2) + Math.Pow((z1 - z2), 2)));
            double cosxX = -(x1 - x2) / Lx;
            double cosxY = -(y1 - y2) / Lx;
            double cosxZ = -(z1 - z2) / Lx;
            double Ly = Math.Sqrt((Math.Pow(((y1 - y2) * ((x1 - x2) * (y1 - y3) - (x1 - x3) * (y1 - y2)) + (z1 - z2) * ((x1 - x2) * (z1 - z3) - (x1 - x3) * (z1 - z2))), 2) + Math.Pow(((x1 - x2) * ((x1 - x2) * (y1 - y3) - (x1 - x3) * (y1 - y2)) - (z1 - z2) * ((y1 - y2) * (z1 - z3) - (y1 - y3) * (z1 - z2))), 2) + Math.Pow(((x1 - x2) * ((x1 - x2) * (z1 - z3) - (x1 - x3) * (z1 - z2)) + (y1 - y2) * ((y1 - y2) * (z1 - z3) - (y1 - y3) * (z1 - z2))), 2)));
            double cosyX = ((y1 - y2) * ((x1 - x2) * (y1 - y3) - (x1 - x3) * (y1 - y2)) + (z1 - z2) * ((x1 - x2) * (z1 - z3) - (x1 - x3) * (z1 - z2))) / Ly;
            double cosyY = -((x1 - x2) * ((x1 - x2) * (y1 - y3) - (x1 - x3) * (y1 - y2)) - (z1 - z2) * ((y1 - y2) * (z1 - z3) - (y1 - y3) * (z1 - z2))) / Ly;
            double cosyZ = -((x1 - x2) * ((x1 - x2) * (z1 - z3) - (x1 - x3) * (z1 - z2)) + (y1 - y2) * ((y1 - y2) * (z1 - z3) - (y1 - y3) * (z1 - z2))) / Ly;
            double Lz = Math.Sqrt((Math.Pow(((x1 - x2) * (y1 - y3) - (x1 - x3) * (y1 - y2)), 2) + Math.Pow(((x1 - x2) * (z1 - z3) - (x1 - x3) * (z1 - z2)), 2) + Math.Pow(((y1 - y2) * (z1 - z3) - (y1 - y3) * (z1 - z2)), 2)));
            double coszX = ((y1 - y2) * (z1 - z3) - (y1 - y3) * (z1 - z2)) / Lz;
            double coszY = -((x1 - x2) * (z1 - z3) - (x1 - x3) * (z1 - z2)) / Lz;
            double coszZ = ((x1 - x2) * (y1 - y3) - (x1 - x3) * (y1 - y2)) / Lz;

            // assembling nodal x,y,z tranformation matrix tf
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

            // assemble the full transformation matrix T for the entire element (12x12 matrix)
            Matrix<double> one = Matrix<double>.Build.Dense(1, 1, 1);
            var T = tf;
            T = T.DiagonalStack(one);
            T = T.DiagonalStack(tf);
            T = T.DiagonalStack(one);
            T = T.DiagonalStack(tf);
            T = T.DiagonalStack(one);
            Matrix<double> T_T = T.Transpose(); // and the transposed tranformation matrix

            // initiates the local coordinate matrix, initiated with global coordinates
            Matrix<double> lcoord = Matrix<double>.Build.DenseOfArray(new double[,]
            {
                { x1, x2, x3 },
                { y1, y2, y3 },
                { z1, z2, z3 }
            });

            //transforms lcoord into local coordinate values
            lcoord = tf.Multiply(lcoord);

            // sets the new (local) coordinate values
            x1 = lcoord[0, 0];
            x2 = lcoord[0, 1];
            x3 = lcoord[0, 2];
            y1 = lcoord[1, 0];
            y2 = lcoord[1, 1];
            y3 = lcoord[1, 2];
            z1 = lcoord[2, 0];
            z2 = lcoord[2, 1];
            z3 = lcoord[2, 2]; // Note that z1 = z2 = z3, if all goes according to plan

            #endregion

            double Area = Math.Abs(0.5 * (x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2)));

            // Establishes the general flexural rigidity matrix for plate
            Matrix<double> C = Matrix<double>.Build.Dense(3, 3);
            C[0, 0] = 1;
            C[0, 1] = nu;
            C[1, 0] = nu;
            C[1, 1] = 1;
            C[2, 2] = (1 - nu)*0.5;

            double C_add = E / (1 - Math.Pow(nu, 2)); // additional part to add to every indice in C matrix

            #region Morley Bending Triangle -- Bending part of element gives [z1 z2 z3 phi1 phi2 phi3]

            Matrix<double> lcoord_temp = Matrix<double>.Build.DenseOfArray(new double[,] { { x1 }, { y1 }, { z1 } });
            lcoord = lcoord.Append(lcoord_temp);

            //// defines variables for simplicity
            double x13 = x1 - x3;
            double x32 = x3 - x2;
            double y23 = y2 - y3;
            double y31 = y3 - y1;

            double[] ga = new double[3];
            double[] my = new double[3];
            double[] a = new double[3];

            for (int i = 0; i < 3; i++)
            {
                double c, s;
                double len = Math.Sqrt(Math.Pow(lcoord[0, i + 1] - lcoord[0, i], 2) + Math.Pow(lcoord[1, i + 1] - lcoord[1, i], 2));
                if (lcoord[0, i + 1] > lcoord[0, i])
                {
                    c = (lcoord[0, i + 1] - lcoord[0, i]) / len;
                    s = (lcoord[1, i + 1] - lcoord[1, i]) / len;
                }
                else if (lcoord[0, i + 1] < lcoord[0, i])
                {
                    c = (lcoord[0, i] - lcoord[0, i + 1]) / len;
                    s = (lcoord[1, i] - lcoord[1, i + 1]) / len;
                }
                else
                {
                    c = 0.0;
                    s = 1.0;
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

            Matrix<double> Bk_b = Matrix<double>.Build.Dense(3, 6); // Exported from Matlab

            Bk_b[0, 0] = 2 * Math.Pow(y23, 2) - (2 * ga6 * y23 * y31) / my6;
            Bk_b[0, 1] = 2 * Math.Pow(y31, 2) - (2 * my5 * y23 * y31) / ga5;
            Bk_b[0, 2] = 2 * y23 * y31 * (a5 / ga5 + a6 / my6);
            Bk_b[0, 4] = (2 * y23 * y31) / ga5;
            Bk_b[0, 5] = (2 * y23 * y31) / my6;
            Bk_b[1, 0] = 2 * Math.Pow(x32, 2) - (2 * ga6 * x13 * x32) / my6;
            Bk_b[1, 1] = 2 * Math.Pow(x13, 2) - (2 * my5 * x13 * x32) / ga5;
            Bk_b[1, 2] = 2 * x13 * x32 * (a5 / ga5 + a6 / my6);
            Bk_b[1, 4] = (2 * x13 * x32) / ga5;
            Bk_b[1, 5] = (2 * x13 * x32) / my6;
            Bk_b[2, 0] = 4 * x32 * y23 - (ga6 * (2 * x13 * y23 + 2 * x32 * y31)) / my6;
            Bk_b[2, 1] = 4 * x13 * y31 - (my5 * (2 * x13 * y23 + 2 * x32 * y31)) / ga5;
            Bk_b[2, 2] = (2 * x13 * y23 + 2 * x32 * y31) * (a5 / ga5 + a6 / my6);
            Bk_b[2, 4] = (2 * x13 * y23 + 2 * x32 * y31) / ga5;
            Bk_b[2, 5] = (2 * x13 * y23 + 2 * x32 * y31) / my6;

            double Bk_b_add = 1 / (4.0 * Math.Pow(Area, 2)); // additional part to add to every indice in B matrix

            Matrix<double> Bk_b_T = Bk_b.Transpose();

            Matrix<double> ke_b = C.Multiply(Bk_b); // the bending part of the element stiffness matrix
            ke_b = Bk_b_T.Multiply(ke_b);
            double ke_b_add = (Area * t * t * t) / 12; // additional part to add to every indice in ke_b matrix
            ke_b_add = ke_b_add * Bk_b_add * C_add * Bk_b_add; // multiply upp all additional parts
            ke_b = ke_b.Multiply(ke_b_add);

#endregion


            #region Constant Strain/Stress Triangle (CST) -- Membrane part of element gives [x1 y1 x2 y2 x3 y3]

            Matrix<double> Bk_m = Matrix<double>.Build.Dense(3, 6); // Exported from Matlab

            Bk_m[0, 0] = (y2 - y3) / (x1 * y2 - x2 * y1 - x1 * y3 + x3 * y1 + x2 * y3 - x3 * y2);
            Bk_m[0, 2] = -(y1 - y3) / (x1 * y2 - x2 * y1 - x1 * y3 + x3 * y1 + x2 * y3 - x3 * y2);
            Bk_m[0, 4] = (y1 - y2) / (x1 * y2 - x2 * y1 - x1 * y3 + x3 * y1 + x2 * y3 - x3 * y2);
            Bk_m[1, 1] = -(x2 - x3) / (x1 * y2 - x2 * y1 - x1 * y3 + x3 * y1 + x2 * y3 - x3 * y2);
            Bk_m[1, 3] = (x1 - x3) / (x1 * y2 - x2 * y1 - x1 * y3 + x3 * y1 + x2 * y3 - x3 * y2);
            Bk_m[1, 5] = -(x1 - x2) / (x1 * y2 - x2 * y1 - x1 * y3 + x3 * y1 + x2 * y3 - x3 * y2);
            Bk_m[2, 0] = -(x2 - x3) / (x1 * y2 - x2 * y1 - x1 * y3 + x3 * y1 + x2 * y3 - x3 * y2);
            Bk_m[2, 1] = (y2 - y3) / (x1 * y2 - x2 * y1 - x1 * y3 + x3 * y1 + x2 * y3 - x3 * y2);
            Bk_m[2, 2] = (x1 - x3) / (x1 * y2 - x2 * y1 - x1 * y3 + x3 * y1 + x2 * y3 - x3 * y2);
            Bk_m[2, 3] = -(y1 - y3) / (x1 * y2 - x2 * y1 - x1 * y3 + x3 * y1 + x2 * y3 - x3 * y2);
            Bk_m[2, 4] = -(x1 - x2) / (x1 * y2 - x2 * y1 - x1 * y3 + x3 * y1 + x2 * y3 - x3 * y2);
            Bk_m[2, 5] = (y1 - y2) / (x1 * y2 - x2 * y1 - x1 * y3 + x3 * y1 + x2 * y3 - x3 * y2);

            Matrix<double> Bk_m_T = Bk_m.Transpose();

            Matrix<double> ke_m = C.Multiply(Bk_m); // the membrane part of the element stiffness matrix
            ke_m = Bk_m_T.Multiply(ke_m);
            ke_m = ke_m.Multiply(C_add * Area * t);

            #endregion

            // input membrane and bending part into full element stiffness matrix
            // and stacking them from [x1 y1 x2 y2 x3 y3 z1 z2 z3 phi1 phi2 phi3]
            // into [x1 y1 z1 phi1 x2 y2 z2 phi2 x3 y3 z3 phi3] which gives the stacking order: { 0 1 6 9 2 3 7 10 4 5 8 11 }
            Matrix<double> ke = ke_m.DiagonalStack(ke_b);
            ke = SymmetricRearrangeMatrix(ke, new int[] { 0, 1, 6, 9, 2, 3, 7, 10, 4, 5, 8, 11 });


            #region Directly calculate ke
            //ke calculated in matlab script Simplest_shell_triangle.m in local xy coordinates
            //Matrix<double> ke = Matrix<double>.Build.Dense(12, 12);
            //ke[0, 0] = -(Area * E * t * (Math.Pow(x2, 2) - 4 * y2 * y3 - nu * Math.Pow(x2, 2) - nu * Math.Pow(x3, 2) - 2 * x2 * x3 + Math.Pow(x3, 2) + 2 * Math.Pow(y2, 2) + 2 * Math.Pow(y3, 2) + 2 * nu * x2 * x3)) / (2 * (Math.Pow(nu, 2) - 1) * Math.Pow((x1 * y2 - x2 * y1 - x1 * y3 + x3 * y1 + x2 * y3 - x3 * y2), 2));
            //ke[0, 1] = (Area * E * t * (x2 - x3) * (y2 - y3) * (nu + 1)) / (2 * (Math.Pow(nu, 2) - 1) * Math.Pow((x1 * y2 - x2 * y1 - x1 * y3 + x3 * y1 + x2 * y3 - x3 * y2), 2));
            //ke[0, 4] = (Area * E * t * (x1 * x2 - x1 * x3 - x2 * x3 + 2 * y1 * y2 - 2 * y1 * y3 - 2 * y2 * y3 - nu * Math.Pow(x3, 2) + Math.Pow(x3, 2) + 2 * Math.Pow(y3, 2) - nu * x1 * x2 + nu * x1 * x3 + nu * x2 * x3)) / (2 * (Math.Pow(nu, 2) - 1) * Math.Pow((x1 * y2 - x2 * y1 - x1 * y3 + x3 * y1 + x2 * y3 - x3 * y2), 2));
            //ke[0, 5] = -(Area * E * t * (x2 * y1 - x3 * y1 - x2 * y3 + x3 * y3 + 2 * nu * x1 * y2 - nu * x2 * y1 - 2 * nu * x1 * y3 + nu * x3 * y1 + nu * x2 * y3 - 2 * nu * x3 * y2 + nu * x3 * y3)) / (2 * (Math.Pow(nu, 2) - 1) * Math.Pow((x1 * y2 - x2 * y1 - x1 * y3 + x3 * y1 + x2 * y3 - x3 * y2), 2));
            //ke[0, 8] = -(Area * E * t * (x1 * x2 - x1 * x3 + x2 * x3 + 2 * y1 * y2 - 2 * y1 * y3 + 2 * y2 * y3 + nu * Math.Pow(x2, 2) - Math.Pow(x2, 2) - 2 * Math.Pow(y2, 2) - nu * x1 * x2 + nu * x1 * x3 - nu * x2 * x3)) / (2 * (Math.Pow(nu, 2) - 1) * Math.Pow((x1 * y2 - x2 * y1 - x1 * y3 + x3 * y1 + x2 * y3 - x3 * y2), 2));
            //ke[0, 9] = -(Area * E * t * (x2 * y2 - x2 * y1 + x3 * y1 - x3 * y2 - 2 * nu * x1 * y2 + nu * x2 * y1 + 2 * nu * x1 * y3 + nu * x2 * y2 - nu * x3 * y1 - 2 * nu * x2 * y3 + nu * x3 * y2)) / (2 * (Math.Pow(nu, 2) - 1) * Math.Pow((x1 * y2 - x2 * y1 - x1 * y3 + x3 * y1 + x2 * y3 - x3 * y2), 2));
            //ke[1, 0] = (Area * E * t * (x2 - x3) * (y2 - y3) * (nu + 1)) / (2 * (Math.Pow(nu, 2) - 1) * Math.Pow((x1 * y2 - x2 * y1 - x1 * y3 + x3 * y1 + x2 * y3 - x3 * y2), 2));
            //ke[1, 1] = -(Area * E * t * (2 * Math.Pow(x2, 2) - 2 * y2 * y3 - nu * Math.Pow(y2, 2) - nu * Math.Pow(y3, 2) - 4 * x2 * x3 + 2 * Math.Pow(x3, 2) + Math.Pow(y2, 2) + Math.Pow(y3, 2) + 2 * nu * y2 * y3)) / (2 * (Math.Pow(nu, 2) - 1) * Math.Pow((x1 * y2 - x2 * y1 - x1 * y3 + x3 * y1 + x2 * y3 - x3 * y2), 2));
            //ke[1, 4] = -(Area * E * t * (x1 * y2 - x1 * y3 - x3 * y2 + x3 * y3 - nu * x1 * y2 + 2 * nu * x2 * y1 + nu * x1 * y3 - 2 * nu * x3 * y1 - 2 * nu * x2 * y3 + nu * x3 * y2 + nu * x3 * y3)) / (2 * (Math.Pow(nu, 2) - 1) * Math.Pow((x1 * y2 - x2 * y1 - x1 * y3 + x3 * y1 + x2 * y3 - x3 * y2), 2));
            //ke[1, 5] = (Area * E * t * (2 * x1 * x2 - 2 * x1 * x3 - 2 * x2 * x3 + y1 * y2 - y1 * y3 - y2 * y3 - nu * Math.Pow(y3, 2) + 2 * Math.Pow(x3, 2) + Math.Pow(y3, 2) - nu * y1 * y2 + nu * y1 * y3 + nu * y2 * y3)) / (2 * (Math.Pow(nu, 2) - 1) * Math.Pow((x1 * y2 - x2 * y1 - x1 * y3 + x3 * y1 + x2 * y3 - x3 * y2), 2));
            //ke[1, 8] = -(Area * E * t * (x1 * y3 - x1 * y2 + x2 * y2 - x2 * y3 + nu * x1 * y2 - 2 * nu * x2 * y1 - nu * x1 * y3 + nu * x2 * y2 + 2 * nu * x3 * y1 + nu * x2 * y3 - 2 * nu * x3 * y2)) / (2 * (Math.Pow(nu, 2) - 1) * Math.Pow((x1 * y2 - x2 * y1 - x1 * y3 + x3 * y1 + x2 * y3 - x3 * y2), 2));
            //ke[1, 9] = -(Area * E * t * (2 * x1 * x2 - 2 * x1 * x3 + 2 * x2 * x3 + y1 * y2 - y1 * y3 + y2 * y3 + nu * Math.Pow(y2, 2) - 2 * Math.Pow(x2, 2) - Math.Pow(y2, 2) - nu * y1 * y2 + nu * y1 * y3 - nu * y2 * y3)) / (2 * (Math.Pow(nu, 2) - 1) * Math.Pow((x1 * y2 - x2 * y1 - x1 * y3 + x3 * y1 + x2 * y3 - x3 * y2), 2));
            //ke[2, 2] = -(E * Math.Pow(t, 3) * ((2 * Math.Pow(x32, 2) - (2 * ga6 * x13 * x32) / my6) * (2 * Math.Pow(x32, 2) + nu * (2 * Math.Pow(y23, 2) - (2 * ga6 * y23 * y31) / my6) - (2 * ga6 * x13 * x32) / my6) - (nu / 2 - 1 / 2) * Math.Pow((4 * x32 * y23 - (ga6 * (2 * x13 * y23 + 2 * x32 * y31)) / my6), 2) + (2 * Math.Pow(y23, 2) - (2 * ga6 * y23 * y31) / my6) * (2 * Math.Pow(y23, 2) + nu * (2 * Math.Pow(x32, 2) - (2 * ga6 * x13 * x32) / my6) - (2 * ga6 * y23 * y31) / my6))) / (192 * Math.Pow(Area, 3) * (Math.Pow(nu, 2) - 1));
            //ke[2, 6] = -(E * Math.Pow(t, 3) * ((2 * Math.Pow(x13, 2) - (2 * my5 * x13 * x32) / ga5) * (2 * Math.Pow(x32, 2) + nu * (2 * Math.Pow(y23, 2) - (2 * ga6 * y23 * y31) / my6) - (2 * ga6 * x13 * x32) / my6) + (2 * Math.Pow(y31, 2) - (2 * my5 * y23 * y31) / ga5) * (2 * Math.Pow(y23, 2) + nu * (2 * Math.Pow(x32, 2) - (2 * ga6 * x13 * x32) / my6) - (2 * ga6 * y23 * y31) / my6) - (nu / 2 - 1 / 2) * (4 * x13 * y31 - (my5 * (2 * x13 * y23 + 2 * x32 * y31)) / ga5) * (4 * x32 * y23 - (ga6 * (2 * x13 * y23 + 2 * x32 * y31)) / my6))) / (192 * Math.Pow(Area, 3) * (Math.Pow(nu, 2) - 1));
            //ke[2, 7] = -(E * Math.Pow(t, 3) * ((2 * x13 * x32 * (2 * Math.Pow(x32, 2) + nu * (2 * Math.Pow(y23, 2) - (2 * ga6 * y23 * y31) / my6) - (2 * ga6 * x13 * x32) / my6)) / ga5 + (2 * y23 * y31 * (2 * Math.Pow(y23, 2) + nu * (2 * Math.Pow(x32, 2) - (2 * ga6 * x13 * x32) / my6) - (2 * ga6 * y23 * y31) / my6)) / ga5 - ((2 * x13 * y23 + 2 * x32 * y31) * (nu / 2 - 1 / 2) * (4 * x32 * y23 - (ga6 * (2 * x13 * y23 + 2 * x32 * y31)) / my6)) / ga5)) / (192 * Math.Pow(Area, 3) * (Math.Pow(nu, 2) - 1));
            //ke[2, 10] = -(E * Math.Pow(t, 3) * (2 * x13 * x32 * (a5 / ga5 + a6 / my6) * (2 * Math.Pow(x32, 2) + nu * (2 * Math.Pow(y23, 2) - (2 * ga6 * y23 * y31) / my6) - (2 * ga6 * x13 * x32) / my6) - (2 * x13 * y23 + 2 * x32 * y31) * (nu / 2 - 1 / 2) * (4 * x32 * y23 - (ga6 * (2 * x13 * y23 + 2 * x32 * y31)) / my6) * (a5 / ga5 + a6 / my6) + 2 * y23 * y31 * (a5 / ga5 + a6 / my6) * (2 * Math.Pow(y23, 2) + nu * (2 * Math.Pow(x32, 2) - (2 * ga6 * x13 * x32) / my6) - (2 * ga6 * y23 * y31) / my6))) / (192 * Math.Pow(Area, 3) * (Math.Pow(nu, 2) - 1));
            //ke[2, 11] = -(E * Math.Pow(t, 3) * ((4 * x13 * x32 * (my6 * Math.Pow(x32, 2) + my6 * nu * Math.Pow(y23, 2) - ga6 * x13 * x32 - ga6 * nu * y23 * y31)) / Math.Pow(my6, 2) + (4 * y23 * y31 * (my6 * Math.Pow(y23, 2) + my6 * nu * Math.Pow(x32, 2) - ga6 * y23 * y31 - ga6 * nu * x13 * x32)) / Math.Pow(my6, 2) + (2 * (2 * x13 * y23 + 2 * x32 * y31) * (nu / 2 - 1 / 2) * (ga6 * x13 * y23 + ga6 * x32 * y31 - 2 * my6 * x32 * y23)) / Math.Pow(my6, 2))) / (192 * Math.Pow(Area, 3) * (Math.Pow(nu, 2) - 1));
            //ke[4, 0] = (Area * E * t * (x1 * x2 - x1 * x3 - x2 * x3 + 2 * y1 * y2 - 2 * y1 * y3 - 2 * y2 * y3 - nu * Math.Pow(x3, 2) + Math.Pow(x3, 2) + 2 * Math.Pow(y3, 2) - nu * x1 * x2 + nu * x1 * x3 + nu * x2 * x3)) / (2 * (Math.Pow(nu, 2) - 1) * Math.Pow((x1 * y2 - x2 * y1 - x1 * y3 + x3 * y1 + x2 * y3 - x3 * y2), 2));
            //ke[4, 1] = -(Area * E * t * (x1 * y2 - x1 * y3 - x3 * y2 + x3 * y3 - nu * x1 * y2 + 2 * nu * x2 * y1 + nu * x1 * y3 - 2 * nu * x3 * y1 - 2 * nu * x2 * y3 + nu * x3 * y2 + nu * x3 * y3)) / (2 * (Math.Pow(nu, 2) - 1) * Math.Pow((x1 * y2 - x2 * y1 - x1 * y3 + x3 * y1 + x2 * y3 - x3 * y2), 2));
            //ke[4, 4] = -(Area * E * t * (Math.Pow(x1, 2) - 4 * y1 * y3 - nu * Math.Pow(x1, 2) - nu * Math.Pow(x3, 2) - 2 * x1 * x3 + Math.Pow(x3, 2) + 2 * Math.Pow(y1, 2) + 2 * Math.Pow(y3, 2) + 2 * nu * x1 * x3)) / (2 * (Math.Pow(nu, 2) - 1) * Math.Pow((x1 * y2 - x2 * y1 - x1 * y3 + x3 * y1 + x2 * y3 - x3 * y2), 2));
            //ke[4, 5] = (Area * E * t * (x1 - x3) * (y1 - y3) * (nu + 1)) / (2 * (Math.Pow(nu, 2) - 1) * Math.Pow((x1 * y2 - x2 * y1 - x1 * y3 + x3 * y1 + x2 * y3 - x3 * y2), 2));
            //ke[4, 8] = -(Area * E * t * (x1 * x2 + x1 * x3 - x2 * x3 + 2 * y1 * y2 + 2 * y1 * y3 - 2 * y2 * y3 + nu * Math.Pow(x1, 2) - Math.Pow(x1, 2) - 2 * Math.Pow(y1, 2) - nu * x1 * x2 - nu * x1 * x3 + nu * x2 * x3)) / (2 * (Math.Pow(nu, 2) - 1) * Math.Pow((x1 * y2 - x2 * y1 - x1 * y3 + x3 * y1 + x2 * y3 - x3 * y2), 2));
            //ke[4, 9] = -(Area * E * t * (x1 * y1 - x1 * y2 - x3 * y1 + x3 * y2 + nu * x1 * y1 + nu * x1 * y2 - 2 * nu * x2 * y1 - 2 * nu * x1 * y3 + nu * x3 * y1 + 2 * nu * x2 * y3 - nu * x3 * y2)) / (2 * (Math.Pow(nu, 2) - 1) * Math.Pow((x1 * y2 - x2 * y1 - x1 * y3 + x3 * y1 + x2 * y3 - x3 * y2), 2));
            //ke[5, 0] = -(Area * E * t * (x2 * y1 - x3 * y1 - x2 * y3 + x3 * y3 + 2 * nu * x1 * y2 - nu * x2 * y1 - 2 * nu * x1 * y3 + nu * x3 * y1 + nu * x2 * y3 - 2 * nu * x3 * y2 + nu * x3 * y3)) / (2 * (Math.Pow(nu, 2) - 1) * Math.Pow((x1 * y2 - x2 * y1 - x1 * y3 + x3 * y1 + x2 * y3 - x3 * y2), 2));
            //ke[5, 1] = (Area * E * t * (2 * x1 * x2 - 2 * x1 * x3 - 2 * x2 * x3 + y1 * y2 - y1 * y3 - y2 * y3 - nu * Math.Pow(y3, 2) + 2 * Math.Pow(x3, 2) + Math.Pow(y3, 2) - nu * y1 * y2 + nu * y1 * y3 + nu * y2 * y3)) / (2 * (Math.Pow(nu, 2) - 1) * Math.Pow((x1 * y2 - x2 * y1 - x1 * y3 + x3 * y1 + x2 * y3 - x3 * y2), 2));
            //ke[5, 4] = (Area * E * t * (x1 - x3) * (y1 - y3) * (nu + 1)) / (2 * (Math.Pow(nu, 2) - 1) * Math.Pow((x1 * y2 - x2 * y1 - x1 * y3 + x3 * y1 + x2 * y3 - x3 * y2), 2));
            //ke[5, 5] = -(Area * E * t * (2 * Math.Pow(x1, 2) - 2 * y1 * y3 - nu * Math.Pow(y1, 2) - nu * Math.Pow(y3, 2) - 4 * x1 * x3 + 2 * Math.Pow(x3, 2) + Math.Pow(y1, 2) + Math.Pow(y3, 2) + 2 * nu * y1 * y3)) / (2 * (Math.Pow(nu, 2) - 1) * Math.Pow((x1 * y2 - x2 * y1 - x1 * y3 + x3 * y1 + x2 * y3 - x3 * y2), 2));
            //ke[5, 8] = -(Area * E * t * (x1 * y1 - x2 * y1 - x1 * y3 + x2 * y3 + nu * x1 * y1 - 2 * nu * x1 * y2 + nu * x2 * y1 + nu * x1 * y3 - 2 * nu * x3 * y1 - nu * x2 * y3 + 2 * nu * x3 * y2)) / (2 * (Math.Pow(nu, 2) - 1) * Math.Pow((x1 * y2 - x2 * y1 - x1 * y3 + x3 * y1 + x2 * y3 - x3 * y2), 2));
            //ke[5, 9] = -(Area * E * t * (2 * x1 * x2 + 2 * x1 * x3 - 2 * x2 * x3 + y1 * y2 + y1 * y3 - y2 * y3 + nu * Math.Pow(y1, 2) - 2 * Math.Pow(x1, 2) - Math.Pow(y1, 2) - nu * y1 * y2 - nu * y1 * y3 + nu * y2 * y3)) / (2 * (Math.Pow(nu, 2) - 1) * Math.Pow((x1 * y2 - x2 * y1 - x1 * y3 + x3 * y1 + x2 * y3 - x3 * y2), 2));
            //ke[6, 2] = -(E * Math.Pow(t, 3) * ((2 * Math.Pow(x32, 2) - (2 * ga6 * x13 * x32) / my6) * (2 * Math.Pow(x13, 2) + nu * (2 * Math.Pow(y31, 2) - (2 * my5 * y23 * y31) / ga5) - (2 * my5 * x13 * x32) / ga5) + (2 * Math.Pow(y23, 2) - (2 * ga6 * y23 * y31) / my6) * (2 * Math.Pow(y31, 2) + nu * (2 * Math.Pow(x13, 2) - (2 * my5 * x13 * x32) / ga5) - (2 * my5 * y23 * y31) / ga5) - (nu / 2 - 1 / 2) * (4 * x13 * y31 - (my5 * (2 * x13 * y23 + 2 * x32 * y31)) / ga5) * (4 * x32 * y23 - (ga6 * (2 * x13 * y23 + 2 * x32 * y31)) / my6))) / (192 * Math.Pow(Area, 3) * (Math.Pow(nu, 2) - 1));
            //ke[6, 6] = -(E * Math.Pow(t, 3) * ((2 * Math.Pow(x13, 2) - (2 * my5 * x13 * x32) / ga5) * (2 * Math.Pow(x13, 2) + nu * (2 * Math.Pow(y31, 2) - (2 * my5 * y23 * y31) / ga5) - (2 * my5 * x13 * x32) / ga5) - (nu / 2 - 1 / 2) * Math.Pow((4 * x13 * y31 - (my5 * (2 * x13 * y23 + 2 * x32 * y31)) / ga5), 2) + (2 * Math.Pow(y31, 2) - (2 * my5 * y23 * y31) / ga5) * (2 * Math.Pow(y31, 2) + nu * (2 * Math.Pow(x13, 2) - (2 * my5 * x13 * x32) / ga5) - (2 * my5 * y23 * y31) / ga5))) / (192 * Math.Pow(Area, 3) * (Math.Pow(nu, 2) - 1));
            //ke[6, 7] = -(E * Math.Pow(t, 3) * ((4 * x13 * x32 * (ga5 * Math.Pow(x13, 2) + ga5 * nu * Math.Pow(y31, 2) - my5 * x13 * x32 - my5 * nu * y23 * y31)) / Math.Pow(ga5, 2) + (4 * y23 * y31 * (ga5 * Math.Pow(y31, 2) + ga5 * nu * Math.Pow(x13, 2) - my5 * y23 * y31 - my5 * nu * x13 * x32)) / Math.Pow(ga5, 2) + (2 * (2 * x13 * y23 + 2 * x32 * y31) * (nu / 2 - 1 / 2) * (my5 * x13 * y23 - 2 * ga5 * x13 * y31 + my5 * x32 * y31)) / Math.Pow(ga5, 2))) / (192 * Math.Pow(Area, 3) * (Math.Pow(nu, 2) - 1));
            //ke[6, 10] = -(E * Math.Pow(t, 3) * (2 * x13 * x32 * (a5 / ga5 + a6 / my6) * (2 * Math.Pow(x13, 2) + nu * (2 * Math.Pow(y31, 2) - (2 * my5 * y23 * y31) / ga5) - (2 * my5 * x13 * x32) / ga5) - (2 * x13 * y23 + 2 * x32 * y31) * (nu / 2 - 1 / 2) * (4 * x13 * y31 - (my5 * (2 * x13 * y23 + 2 * x32 * y31)) / ga5) * (a5 / ga5 + a6 / my6) + 2 * y23 * y31 * (a5 / ga5 + a6 / my6) * (2 * Math.Pow(y31, 2) + nu * (2 * Math.Pow(x13, 2) - (2 * my5 * x13 * x32) / ga5) - (2 * my5 * y23 * y31) / ga5))) / (192 * Math.Pow(Area, 3) * (Math.Pow(nu, 2) - 1));
            //ke[6, 11] = -(E * Math.Pow(t, 3) * ((2 * x13 * x32 * (2 * Math.Pow(x13, 2) + nu * (2 * Math.Pow(y31, 2) - (2 * my5 * y23 * y31) / ga5) - (2 * my5 * x13 * x32) / ga5)) / my6 + (2 * y23 * y31 * (2 * Math.Pow(y31, 2) + nu * (2 * Math.Pow(x13, 2) - (2 * my5 * x13 * x32) / ga5) - (2 * my5 * y23 * y31) / ga5)) / my6 - ((2 * x13 * y23 + 2 * x32 * y31) * (nu / 2 - 1 / 2) * (4 * x13 * y31 - (my5 * (2 * x13 * y23 + 2 * x32 * y31)) / ga5)) / my6)) / (192 * Math.Pow(Area, 3) * (Math.Pow(nu, 2) - 1));
            //ke[7, 2] = (E * Math.Pow(t, 3) * ((4 * x32 * (x13 * x32 + nu * y23 * y31) * (ga6 * x13 - my6 * x32)) / (ga5 * my6) - (2 * (2 * x13 * y23 + 2 * x32 * y31) * (nu / 2 - 1 / 2) * (ga6 * x13 * y23 + ga6 * x32 * y31 - 2 * my6 * x32 * y23)) / (ga5 * my6) + (4 * y23 * (y23 * y31 + nu * x13 * x32) * (ga6 * y31 - my6 * y23)) / (ga5 * my6))) / (192 * Math.Pow(Area, 3) * (Math.Pow(nu, 2) - 1));
            //ke[7, 6] = -(E * Math.Pow(t, 3) * ((2 * (2 * x13 * y23 + 2 * x32 * y31) * (nu / 2 - 1 / 2) * (my5 * x13 * y23 - 2 * ga5 * x13 * y31 + my5 * x32 * y31)) / Math.Pow(ga5, 2) + (4 * x13 * (x13 * x32 + nu * y23 * y31) * (ga5 * x13 - my5 * x32)) / Math.Pow(ga5, 2) + (4 * y31 * (y23 * y31 + nu * x13 * x32) * (ga5 * y31 - my5 * y23)) / Math.Pow(ga5, 2))) / (192 * Math.Pow(Area, 3) * (Math.Pow(nu, 2) - 1));
            //ke[7, 7] = -(E * Math.Pow(t, 3) * ((2 * x13 * x32 * ((2 * x13 * x32) / ga5 + (2 * nu * y23 * y31) / ga5)) / ga5 - (Math.Pow((2 * x13 * y23 + 2 * x32 * y31), 2) * (nu / 2 - 1 / 2)) / Math.Pow(ga5, 2) + (2 * y23 * y31 * ((2 * y23 * y31) / ga5 + (2 * nu * x13 * x32) / ga5)) / ga5)) / (192 * Math.Pow(Area, 3) * (Math.Pow(nu, 2) - 1));
            //ke[7, 10] = -(E * Math.Pow(t, 3) * (a6 * ga5 + a5 * my6) * (2 * Math.Pow(x13, 2) * Math.Pow(x32, 2) + Math.Pow(x13, 2) * Math.Pow(y23, 2) + Math.Pow(x32, 2) * Math.Pow(y31, 2) + 2 * Math.Pow(y23, 2) * Math.Pow(y31, 2) - nu * Math.Pow(x13, 2) * Math.Pow(y23, 2) - nu * Math.Pow(x32, 2) * Math.Pow(y31, 2) + 2 * x13 * x32 * y23 * y31 + 2 * nu * x13 * x32 * y23 * y31)) / (96 * Math.Pow(Area, 3) * Math.Pow(ga5, 2) * my6 * (Math.Pow(nu, 2) - 1));
            //ke[7, 11] = -(E * Math.Pow(t, 3) * (2 * Math.Pow(x13, 2) * Math.Pow(x32, 2) + Math.Pow(x13, 2) * Math.Pow(y23, 2) + Math.Pow(x32, 2) * Math.Pow(y31, 2) + 2 * Math.Pow(y23, 2) * Math.Pow(y31, 2) - nu * Math.Pow(x13, 2) * Math.Pow(y23, 2) - nu * Math.Pow(x32, 2) * Math.Pow(y31, 2) + 2 * x13 * x32 * y23 * y31 + 2 * nu * x13 * x32 * y23 * y31)) / (96 * Math.Pow(Area, 3) * ga5 * my6 * (Math.Pow(nu, 2) - 1));
            //ke[8, 0] = -(Area * E * t * (x1 * x2 - x1 * x3 + x2 * x3 + 2 * y1 * y2 - 2 * y1 * y3 + 2 * y2 * y3 + nu * Math.Pow(x2, 2) - Math.Pow(x2, 2) - 2 * Math.Pow(y2, 2) - nu * x1 * x2 + nu * x1 * x3 - nu * x2 * x3)) / (2 * (Math.Pow(nu, 2) - 1) * Math.Pow((x1 * y2 - x2 * y1 - x1 * y3 + x3 * y1 + x2 * y3 - x3 * y2), 2));
            //ke[8, 1] = -(Area * E * t * (x1 * y3 - x1 * y2 + x2 * y2 - x2 * y3 + nu * x1 * y2 - 2 * nu * x2 * y1 - nu * x1 * y3 + nu * x2 * y2 + 2 * nu * x3 * y1 + nu * x2 * y3 - 2 * nu * x3 * y2)) / (2 * (Math.Pow(nu, 2) - 1) * Math.Pow((x1 * y2 - x2 * y1 - x1 * y3 + x3 * y1 + x2 * y3 - x3 * y2), 2));
            //ke[8, 4] = -(Area * E * t * (x1 * x2 + x1 * x3 - x2 * x3 + 2 * y1 * y2 + 2 * y1 * y3 - 2 * y2 * y3 + nu * Math.Pow(x1, 2) - Math.Pow(x1, 2) - 2 * Math.Pow(y1, 2) - nu * x1 * x2 - nu * x1 * x3 + nu * x2 * x3)) / (2 * (Math.Pow(nu, 2) - 1) * Math.Pow((x1 * y2 - x2 * y1 - x1 * y3 + x3 * y1 + x2 * y3 - x3 * y2), 2));
            //ke[8, 5] = -(Area * E * t * (x1 * y1 - x2 * y1 - x1 * y3 + x2 * y3 + nu * x1 * y1 - 2 * nu * x1 * y2 + nu * x2 * y1 + nu * x1 * y3 - 2 * nu * x3 * y1 - nu * x2 * y3 + 2 * nu * x3 * y2)) / (2 * (Math.Pow(nu, 2) - 1) * Math.Pow((x1 * y2 - x2 * y1 - x1 * y3 + x3 * y1 + x2 * y3 - x3 * y2), 2));
            //ke[8, 8] = -(Area * E * t * (Math.Pow(x1, 2) - 4 * y1 * y2 - nu * Math.Pow(x1, 2) - nu * Math.Pow(x2, 2) - 2 * x1 * x2 + Math.Pow(x2, 2) + 2 * Math.Pow(y1, 2) + 2 * Math.Pow(y2, 2) + 2 * nu * x1 * x2)) / (2 * (Math.Pow(nu, 2) - 1) * Math.Pow((x1 * y2 - x2 * y1 - x1 * y3 + x3 * y1 + x2 * y3 - x3 * y2), 2));
            //ke[8, 9] = (Area * E * t * (x1 - x2) * (y1 - y2) * (nu + 1)) / (2 * (Math.Pow(nu, 2) - 1) * Math.Pow((x1 * y2 - x2 * y1 - x1 * y3 + x3 * y1 + x2 * y3 - x3 * y2), 2));
            //ke[9, 0] = -(Area * E * t * (x2 * y2 - x2 * y1 + x3 * y1 - x3 * y2 - 2 * nu * x1 * y2 + nu * x2 * y1 + 2 * nu * x1 * y3 + nu * x2 * y2 - nu * x3 * y1 - 2 * nu * x2 * y3 + nu * x3 * y2)) / (2 * (Math.Pow(nu, 2) - 1) * Math.Pow((x1 * y2 - x2 * y1 - x1 * y3 + x3 * y1 + x2 * y3 - x3 * y2), 2));
            //ke[9, 1] = -(Area * E * t * (2 * x1 * x2 - 2 * x1 * x3 + 2 * x2 * x3 + y1 * y2 - y1 * y3 + y2 * y3 + nu * Math.Pow(y2, 2) - 2 * Math.Pow(x2, 2) - Math.Pow(y2, 2) - nu * y1 * y2 + nu * y1 * y3 - nu * y2 * y3)) / (2 * (Math.Pow(nu, 2) - 1) * Math.Pow((x1 * y2 - x2 * y1 - x1 * y3 + x3 * y1 + x2 * y3 - x3 * y2), 2));
            //ke[9, 4] = -(Area * E * t * (x1 * y1 - x1 * y2 - x3 * y1 + x3 * y2 + nu * x1 * y1 + nu * x1 * y2 - 2 * nu * x2 * y1 - 2 * nu * x1 * y3 + nu * x3 * y1 + 2 * nu * x2 * y3 - nu * x3 * y2)) / (2 * (Math.Pow(nu, 2) - 1) * Math.Pow((x1 * y2 - x2 * y1 - x1 * y3 + x3 * y1 + x2 * y3 - x3 * y2), 2));
            //ke[9, 5] = -(Area * E * t * (2 * x1 * x2 + 2 * x1 * x3 - 2 * x2 * x3 + y1 * y2 + y1 * y3 - y2 * y3 + nu * Math.Pow(y1, 2) - 2 * Math.Pow(x1, 2) - Math.Pow(y1, 2) - nu * y1 * y2 - nu * y1 * y3 + nu * y2 * y3)) / (2 * (Math.Pow(nu, 2) - 1) * Math.Pow((x1 * y2 - x2 * y1 - x1 * y3 + x3 * y1 + x2 * y3 - x3 * y2), 2));
            //ke[9, 8] = (Area * E * t * (x1 - x2) * (y1 - y2) * (nu + 1)) / (2 * (Math.Pow(nu, 2) - 1) * Math.Pow((x1 * y2 - x2 * y1 - x1 * y3 + x3 * y1 + x2 * y3 - x3 * y2), 2));
            //ke[9, 9] = -(Area * E * t * (2 * Math.Pow(x1, 2) - 2 * y1 * y2 - nu * Math.Pow(y1, 2) - nu * Math.Pow(y2, 2) - 4 * x1 * x2 + 2 * Math.Pow(x2, 2) + Math.Pow(y1, 2) + Math.Pow(y2, 2) + 2 * nu * y1 * y2)) / (2 * (Math.Pow(nu, 2) - 1) * Math.Pow((x1 * y2 - x2 * y1 - x1 * y3 + x3 * y1 + x2 * y3 - x3 * y2), 2));
            //ke[10, 2] = -(E * Math.Pow(t, 3) * ((2 * Math.Pow(x32, 2) - (2 * ga6 * x13 * x32) / my6) * (2 * x13 * x32 * (a5 / ga5 + a6 / my6) + 2 * nu * y23 * y31 * (a5 / ga5 + a6 / my6)) + (2 * Math.Pow(y23, 2) - (2 * ga6 * y23 * y31) / my6) * (2 * y23 * y31 * (a5 / ga5 + a6 / my6) + 2 * nu * x13 * x32 * (a5 / ga5 + a6 / my6)) - (2 * x13 * y23 + 2 * x32 * y31) * (nu / 2 - 1 / 2) * (4 * x32 * y23 - (ga6 * (2 * x13 * y23 + 2 * x32 * y31)) / my6) * (a5 / ga5 + a6 / my6))) / (192 * Math.Pow(Area, 3) * (Math.Pow(nu, 2) - 1));
            //ke[10, 6] = -(E * Math.Pow(t, 3) * ((2 * Math.Pow(x13, 2) - (2 * my5 * x13 * x32) / ga5) * (2 * x13 * x32 * (a5 / ga5 + a6 / my6) + 2 * nu * y23 * y31 * (a5 / ga5 + a6 / my6)) + (2 * Math.Pow(y31, 2) - (2 * my5 * y23 * y31) / ga5) * (2 * y23 * y31 * (a5 / ga5 + a6 / my6) + 2 * nu * x13 * x32 * (a5 / ga5 + a6 / my6)) - (2 * x13 * y23 + 2 * x32 * y31) * (nu / 2 - 1 / 2) * (4 * x13 * y31 - (my5 * (2 * x13 * y23 + 2 * x32 * y31)) / ga5) * (a5 / ga5 + a6 / my6))) / (192 * Math.Pow(Area, 3) * (Math.Pow(nu, 2) - 1));
            //ke[10, 7] = -(E * Math.Pow(t, 3) * (a6 * ga5 + a5 * my6) * (2 * Math.Pow(x13, 2) * Math.Pow(x32, 2) + Math.Pow(x13, 2) * Math.Pow(y23, 2) + Math.Pow(x32, 2) * Math.Pow(y31, 2) + 2 * Math.Pow(y23, 2) * Math.Pow(y31, 2) - nu * Math.Pow(x13, 2) * Math.Pow(y23, 2) - nu * Math.Pow(x32, 2) * Math.Pow(y31, 2) + 2 * x13 * x32 * y23 * y31 + 2 * nu * x13 * x32 * y23 * y31)) / (96 * Math.Pow(Area, 3) * Math.Pow(ga5, 2) * my6 * (Math.Pow(nu, 2) - 1));
            //ke[10, 10] = -(E * Math.Pow(t, 3) * Math.Pow((a6 * ga5 + a5 * my6), 2) * (2 * Math.Pow(x13, 2) * Math.Pow(x32, 2) + Math.Pow(x13, 2) * Math.Pow(y23, 2) + Math.Pow(x32, 2) * Math.Pow(y31, 2) + 2 * Math.Pow(y23, 2) * Math.Pow(y31, 2) - nu * Math.Pow(x13, 2) * Math.Pow(y23, 2) - nu * Math.Pow(x32, 2) * Math.Pow(y31, 2) + 2 * x13 * x32 * y23 * y31 + 2 * nu * x13 * x32 * y23 * y31)) / (96 * Math.Pow(Area, 3) * Math.Pow(ga5, 2) * Math.Pow(my6, 2) * (Math.Pow(nu, 2) - 1));
            //ke[10, 11] = -(E * Math.Pow(t, 3) * (a6 * ga5 + a5 * my6) * (2 * Math.Pow(x13, 2) * Math.Pow(x32, 2) + Math.Pow(x13, 2) * Math.Pow(y23, 2) + Math.Pow(x32, 2) * Math.Pow(y31, 2) + 2 * Math.Pow(y23, 2) * Math.Pow(y31, 2) - nu * Math.Pow(x13, 2) * Math.Pow(y23, 2) - nu * Math.Pow(x32, 2) * Math.Pow(y31, 2) + 2 * x13 * x32 * y23 * y31 + 2 * nu * x13 * x32 * y23 * y31)) / (96 * Math.Pow(Area, 3) * ga5 * Math.Pow(my6, 2) * (Math.Pow(nu, 2) - 1));
            //ke[11, 2] = (E * Math.Pow(t, 3) * ((4 * x32 * (x13 * x32 + nu * y23 * y31) * (ga6 * x13 - my6 * x32)) / Math.Pow(my6, 2) - (2 * (2 * x13 * y23 + 2 * x32 * y31) * (nu / 2 - 1 / 2) * (ga6 * x13 * y23 + ga6 * x32 * y31 - 2 * my6 * x32 * y23)) / Math.Pow(my6, 2) + (4 * y23 * (y23 * y31 + nu * x13 * x32) * (ga6 * y31 - my6 * y23)) / Math.Pow(my6, 2))) / (192 * Math.Pow(Area, 3) * (Math.Pow(nu, 2) - 1));
            //ke[11, 6] = -(E * Math.Pow(t, 3) * ((2 * (2 * x13 * y23 + 2 * x32 * y31) * (nu / 2 - 1 / 2) * (my5 * x13 * y23 - 2 * ga5 * x13 * y31 + my5 * x32 * y31)) / (ga5 * my6) + (4 * x13 * (x13 * x32 + nu * y23 * y31) * (ga5 * x13 - my5 * x32)) / (ga5 * my6) + (4 * y31 * (y23 * y31 + nu * x13 * x32) * (ga5 * y31 - my5 * y23)) / (ga5 * my6))) / (192 * Math.Pow(Area, 3) * (Math.Pow(nu, 2) - 1));
            //ke[11, 7] = -(E * Math.Pow(t, 3) * (2 * Math.Pow(x13, 2) * Math.Pow(x32, 2) + Math.Pow(x13, 2) * Math.Pow(y23, 2) + Math.Pow(x32, 2) * Math.Pow(y31, 2) + 2 * Math.Pow(y23, 2) * Math.Pow(y31, 2) - nu * Math.Pow(x13, 2) * Math.Pow(y23, 2) - nu * Math.Pow(x32, 2) * Math.Pow(y31, 2) + 2 * x13 * x32 * y23 * y31 + 2 * nu * x13 * x32 * y23 * y31)) / (96 * Math.Pow(Area, 3) * ga5 * my6 * (Math.Pow(nu, 2) - 1));
            //ke[11, 10] = -(E * Math.Pow(t, 3) * (a6 * ga5 + a5 * my6) * (2 * Math.Pow(x13, 2) * Math.Pow(x32, 2) + Math.Pow(x13, 2) * Math.Pow(y23, 2) + Math.Pow(x32, 2) * Math.Pow(y31, 2) + 2 * Math.Pow(y23, 2) * Math.Pow(y31, 2) - nu * Math.Pow(x13, 2) * Math.Pow(y23, 2) - nu * Math.Pow(x32, 2) * Math.Pow(y31, 2) + 2 * x13 * x32 * y23 * y31 + 2 * nu * x13 * x32 * y23 * y31)) / (96 * Math.Pow(Area, 3) * ga5 * Math.Pow(my6, 2) * (Math.Pow(nu, 2) - 1));
            //ke[11, 11] = -(E * Math.Pow(t, 3) * ((2 * x13 * x32 * ((2 * x13 * x32) / my6 + (2 * nu * y23 * y31) / my6)) / my6 - (Math.Pow((2 * x13 * y23 + 2 * x32 * y31), 2) * (nu / 2 - 1 / 2)) / Math.Pow(my6, 2) + (2 * y23 * y31 * ((2 * y23 * y31) / my6 + (2 * nu * x13 * x32) / my6)) / my6)) / (192 * Math.Pow(Area, 3) * (Math.Pow(nu, 2) - 1));
            #endregion


            Matrix<double> Ke = ke.Multiply(T);
            Ke = T_T.Multiply(Ke);

            return Ke;
        }

        private Matrix<double> SymmetricRearrangeMatrix(Matrix<double> M, int[] arrangement)
        {
            Matrix<double> M_new = Matrix<double>.Build.Dense(12,12);

            for (int i = 0; i < 12; i++)
            {
                for (int j = 0; j < 12; j++)
                {
                    M_new[i, j] = M[arrangement[i], arrangement[j]];
                }
            }
            return M_new;
        }

        private List<double> CreateLoadList(List<string> loadtxt, List<string> momenttxt, List<Point3d> uniqueNodes, List<MeshFace> faces, List<Point3d> vertices, List<Line> edges)
        {
            //initializing loads with list of doubles of size gdofs and entry values = 0
            List<double> loads = new List<double>(new double[uniqueNodes.Count * 3 + edges.Count]);
            List<double> inputLoads = new List<double>();
            List<Point3d> coordlist = new List<Point3d>();

            //parsing point loads
            for (int i = 0; i < loadtxt.Count; i++)
            {
                string coordstr = (loadtxt[i].Split(':')[0]);
                string loadstr = (loadtxt[i].Split(':')[1]);

                string[] coordstr1 = (coordstr.Split(','));
                string[] loadstr1 = (loadstr.Split(','));

                inputLoads.Add(Math.Round(double.Parse(loadstr1[0]), 2));
                inputLoads.Add(Math.Round(double.Parse(loadstr1[1]), 2));
                inputLoads.Add(Math.Round(double.Parse(loadstr1[2]), 2));

                coordlist.Add(new Point3d(Math.Round(double.Parse(coordstr1[0]), 4), Math.Round(double.Parse(coordstr1[1]), 4), Math.Round(double.Parse(coordstr1[2]), 4)));
            }

            //inputting point loads at correct index in loads list
            foreach (Point3d point in coordlist)
            {
                int gNodeIndex = uniqueNodes.IndexOf(point);
                int lNodeIndex = coordlist.IndexOf(point);
                loads[gNodeIndex * 3 + 0] = inputLoads[lNodeIndex * 3 + 0];
                loads[gNodeIndex * 3 + 1] = inputLoads[lNodeIndex * 3 + 1];
                loads[gNodeIndex * 3 + 2] = inputLoads[lNodeIndex * 3 + 2];
            }
            //resetting variables
            inputLoads.Clear();
            coordlist.Clear();

            #region Moment load input
            ////parsing moment loads
            //for (int i = 0; i < momenttxt.Count; i++) if (momenttxt[0] != "")
            //    {
            //        string coordstr = (momenttxt[i].Split(':')[0]);
            //        string loadstr = (momenttxt[i].Split(':')[1]);

            //        string[] coordstr1 = (coordstr.Split(','));
            //        string[] loadstr1 = (loadstr.Split(','));

            //        inputLoads.Add(Math.Round(double.Parse(loadstr1[0]), 2));
            //        inputLoads.Add(Math.Round(double.Parse(loadstr1[1]), 2));
            //        inputLoads.Add(Math.Round(double.Parse(loadstr1[2]), 2));


            //        coordlist.Add(new Point3d(Math.Round(double.Parse(coordstr1[0]), 2), Math.Round(double.Parse(coordstr1[1]), 2), Math.Round(double.Parse(coordstr1[2]), 2)));
            //    }

            ////inputing moment loads at correct index in loads list
            //    foreach (Point3d point in coordlist)
            //{
            //    int gNodeIndex = uniqueNodes.IndexOf(point);
            //    int lNodeIndex = coordlist.IndexOf(point);
            //    loads[gNodeIndex * ldofs + 3] = inputLoads[lNodeIndex * 3 + 0];
            //    loads[gNodeIndex * ldofs + 4] = inputLoads[lNodeIndex * 3 + 1];
            //    loads[gNodeIndex * ldofs + 5] = inputLoads[lNodeIndex * 3 + 2];
            //}
            #endregion
            return loads;
        }

        private Vector<double> CreateBDCList(List<string> bdctxt, List<Point3d> uniqueNodes, List<MeshFace> faces, List<Point3d> vertices, List<Line> edges)
        {
            //initializing bdc_value as vector of size gdofs, and entry values = 1
            Vector<double> bdc_value = Vector.Build.Dense(uniqueNodes.Count * 3 + edges.Count ,1);
            List<int> bdcs = new List<int>();
            List<Point3d> bdc_points = new List<Point3d>(); //Coordinates relating til bdc_value in for (eg. x y z)
            List<int> fixedRotEdges = new List<int>();
            int rows = bdctxt.Count;

            //Parse string input
            int numOfPoints = 0;
            if ((bdctxt[rows-1] != null) && !bdctxt[rows-1].Contains(":"))
            {
                numOfPoints = rows - 1;
                string[] edgestrtemp = bdctxt[rows - 1].Split(',');
                List<string> edgestr = new List<string>();
                edgestr.AddRange(edgestrtemp);
                for (int i = 0; i < edgestr.Count; i++)
                {
                    fixedRotEdges.Add(int.Parse(edgestr[i]));
                }
            }
            else
            {
                numOfPoints = bdctxt.Count-1;
            }
            for (int i = 0; i < numOfPoints; i++)
            {
                //NB! At the moment, the bdc string either looks like
                //0,0,0:0,0,0
                //or like
                //0,0,0:0,0,0:0,..,n where 0,..,n are the indices of mesh faces containing this vertice, and that should be fixed for rotation
                string coordstr = bdctxt[i].Split(':')[0];
                string bdcstr = bdctxt[i].Split(':')[1];

                string[] coordstr1 = (coordstr.Split(','));
                string[] bdcstr1 = (bdcstr.Split(','));

                bdc_points.Add(new Point3d(Math.Round(double.Parse(coordstr1[0]), 4), Math.Round(double.Parse(coordstr1[1]), 4), Math.Round(double.Parse(coordstr1[2]), 4)));

                bdcs.Add(int.Parse(bdcstr1[0]));
                bdcs.Add(int.Parse(bdcstr1[1]));
                bdcs.Add(int.Parse(bdcstr1[2]));
            }

            //Format to correct entries in bdc_value


            foreach (var point in bdc_points)
            {
                int index = bdc_points.IndexOf(point);
                int i = uniqueNodes.IndexOf(point);
                bdc_value[i * 3 + 0] = bdcs[index * 3 + 0];
                bdc_value[i * 3 + 1] = bdcs[index * 3 + 1];
                bdc_value[i * 3 + 2] = bdcs[index * 3 + 2];
            }
            if (numOfPoints == rows-1)
            {
                foreach (var edgeindex in fixedRotEdges)
                {
                    bdc_value[edgeindex+uniqueNodes.Count*3] = 0;
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

                return Properties.Resources.Calc;
            }
        }

        public override Guid ComponentGuid
        {
            get { return new Guid("3a61d696-911f-46cd-a687-ef48a48575b0"); }
        }

        /// Component Visual//
        public class Attributes_Custom : Grasshopper.Kernel.Attributes.GH_ComponentAttributes
        {
            public Attributes_Custom(GH_Component owner) : base(owner) { }
            protected override void Layout()
            {
                base.Layout();

                Rectangle rec0 = GH_Convert.ToRectangle(Bounds);

                rec0.Height += 22;

                Rectangle rec1 = rec0;
                rec1.X = rec0.Left + 1;
                rec1.Y = rec0.Bottom - 22;
                rec1.Width = (rec0.Width) / 3 + 1;
                rec1.Height = 22;
                rec1.Inflate(-2, -2);

                Rectangle rec2 = rec1;
                rec2.X = rec1.Right + 2;

                Bounds = rec0;
                ButtonBounds = rec1;
                ButtonBounds2 = rec2;

            }

            GH_Palette xColor = GH_Palette.Black;
            GH_Palette yColor = GH_Palette.Grey;

            private Rectangle ButtonBounds { get; set; }
            private Rectangle ButtonBounds2 { get; set; }
            private Rectangle ButtonBounds3 { get; set; }

            protected override void Render(GH_Canvas canvas, Graphics graphics, GH_CanvasChannel channel)
            {
                base.Render(canvas, graphics, channel);
                if (channel == GH_CanvasChannel.Objects)
                {
                    GH_Capsule button;
                    if (startCalc == true)
                    {
                        button = GH_Capsule.CreateTextCapsule(ButtonBounds, ButtonBounds, xColor, "Run: On", 3, 0);
                    }
                    else
                    {
                        button = GH_Capsule.CreateTextCapsule(ButtonBounds, ButtonBounds, yColor, "Run: Off", 3, 0);
                    }
                    button.Render(graphics, Selected, false, false);
                    button.Dispose();
                }
                if (channel == GH_CanvasChannel.Objects)
                {
                    GH_Capsule button2 = GH_Capsule.CreateTextCapsule(ButtonBounds2, ButtonBounds2, yColor, "Run Test", 2, 0);
                    button2.Render(graphics, Selected, Owner.Locked, false);
                    button2.Dispose();
                }
            }

            public override GH_ObjectResponse RespondToMouseDown(GH_Canvas sender, GH_CanvasMouseEvent e)
            {
                if (e.Button == MouseButtons.Left)
                {
                    RectangleF rec = ButtonBounds;
                    if (rec.Contains(e.CanvasLocation))
                    {
                        switchColor("Run");
                        if (xColor == GH_Palette.Black) { setStart("Run", true); Owner.ExpireSolution(true); }
                        if (xColor == GH_Palette.Grey) { setStart("Run", false); Owner.ExpireSolution(true); }
                        
                        sender.Refresh();
                        return GH_ObjectResponse.Handled;
                    }
                    rec = ButtonBounds2;
                    if (rec.Contains(e.CanvasLocation))
                    {
                        switchColor("Run Test");
                        if (yColor == GH_Palette.Black) { setStart("Run Test", true); }
                        if (yColor == GH_Palette.Grey) { setStart("Run Test", false); }
                        sender.Refresh();
                        return GH_ObjectResponse.Handled;
                    }
                }
                return base.RespondToMouseDown(sender, e);
            }

            private void switchColor(string button)
            {
                if (button == "Run")
                {
                    if (xColor == GH_Palette.Black) { xColor = GH_Palette.Grey; }
                    else { xColor = GH_Palette.Black; }
                }
                else if (button == "Run Test")
                {
                    if (yColor == GH_Palette.Black) { yColor = GH_Palette.Grey; }
                    else { yColor = GH_Palette.Black; }
                }
            }
        }
    }
}


