using System;
using System.Collections.Generic;

using Grasshopper.Kernel;
using Rhino.Geometry;
using System.Drawing;
using Grasshopper.GUI.Canvas;
using System.Windows.Forms;
using Grasshopper.GUI;
using MathNet.Numerics.LinearAlgebra;

namespace Shell
{
    public class DeformedGeometry : GH_Component
    {
        public DeformedGeometry()
          : base("DeformedGeometry", "DefG",
              "Description",
              "Koala", "Shell")
        {
        }

        //Initialize startcondition and polynomial order
        static bool startDef = true;

        //Method to allow c hanging of variables via GUI (see Component Visual)
        public static void setToggles(string s, bool i)
        {
            if (s == "Run")
            {
                startDef = i;
            }
        }

        public override void CreateAttributes()
        {
            m_attributes = new Attributes_Custom(this);
        }

        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddNumberParameter("Deformation", "Def", "Deformations from ShellCalc", GH_ParamAccess.list);
            pManager.AddNumberParameter("Stresses", "Stress", "Stresses from ShellCalc", GH_ParamAccess.list);
            pManager.AddMeshParameter("Mesh", "M", "Input Geometry (Mesh format)", GH_ParamAccess.item);
            pManager.AddNumberParameter("Scale", "S", "The Scale Factor for Deformation", GH_ParamAccess.item, 10);
            pManager.AddIntegerParameter("stressdirection", "strdir", "", GH_ParamAccess.item, 0);
        }

        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddMeshParameter("Deformed Geometry", "Def.G.", "Deformed Geometry as mesh", GH_ParamAccess.item);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            if (startDef)
            {
                #region Fetch input
                //Expected inputs and outputs
                List<double> def = new List<double>();
                List<double> stresses = new List<double>();
                Mesh mesh = new Mesh();
                double scale = 10;
                List<Line> defGeometry = new List<Line>();
                List<Point3d> defPoints = new List<Point3d>();

                int dimension = 0;

                //Set expected inputs from Indata
                if (!DA.GetDataList(0, def)) return;
                if (!DA.GetDataList(1, stresses)) return;
                if (!DA.GetData(2, ref mesh)) return;
                if (!DA.GetData(3, ref scale)) return;
                if (!DA.GetData(4, ref dimension)) return;
                #endregion

                #region Decompose Mesh and initiate the new deformed mesh defmesh

                List<Point3d> vertices = new List<Point3d>();
                List<MeshFace> faces = new List<MeshFace>();

                foreach (var vertice in mesh.Vertices)
                {
                    vertices.Add(vertice);
                }
                foreach (var face in mesh.Faces)
                {
                    faces.Add(face);
                }

                Mesh defmesh = new Mesh();

                defmesh.Faces.AddFaces(mesh.Faces); // new mesh without vertices

                #endregion

                #region apply deformations to vertices and add them to defmesh

                List<Point3d> new_vertices = new List<Point3d>(); // list of translated vertices
                int i = 0;

                foreach (var p in vertices)
                {
                    new_vertices.Add(new Point3d(p.X + def[i]*scale, p.Y + def[i + 1]*scale, p.Z + def[i + 2]*scale));
                    i += 3;
                }

                defmesh.Vertices.AddVertices(new_vertices);

                Mesh coloredDefMesh = defmesh.DuplicateMesh();
                SetMeshColors(defmesh, stresses, new_vertices, faces, dimension, out coloredDefMesh);

                #endregion

                //#region Scale deformations
                ////List all nodes (every node only once), numbering them according to list index
                //List<Point3d> points = CreatePointList(geometry);

                //int index = 0;
                ////loops through all points and scales x-, y- and z-dir

                ////u(x) = Na, N = shape func, a = nodal values (dof) 
                //foreach (Point3d point in points)
                //{
                //    //fetch global x,y,z placement of point
                //    double x = point.X;
                //    double y = point.Y;
                //    double z = point.Z;

                //    //scales x and z according to input Scale
                //    defPoints.Add(new Point3d(x + scale * def[index], y + scale * def[index + 1], z + scale * def[index + 2]));
                //    index += 6;
                //}
                //#endregion

                //#region Create geometry
                ////creates deformed geometry based on initial geometry placement
                //foreach (Line line in geometry)
                //{
                //    //fetches index of original start and endpoint
                //    int i1 = points.IndexOf(line.From);
                //    int i2 = points.IndexOf(line.To);

                //    //creates new line based on scaled deformation of said points
                //    defGeometry.Add(new Line(defPoints[i1], defPoints[i2]));
                //}
                //#endregion


                //Set output data
                DA.SetData(0, coloredDefMesh);
            }
        }   //End of main program

        private void SetMeshColors(Mesh meshIn, List<double> stresses, List<Point3d> vertices, List<MeshFace> faces, int direction, out Mesh meshOut)
        {
            meshOut = meshIn.DuplicateMesh();

            double max = 0;
            double min = 0;
            List<int> R = new List<int>(faces.Count);
            List<int> G = new List<int>(faces.Count);
            List<int> B = new List<int>(faces.Count);
            int[,] facesConnectedToVertex = new int[faces.Count,3];

            for (int i = 0; i < stresses.Count/6; i++)
            {
                double stress = stresses[i * 6 + direction];
                if (stress > max)
                {
                    max = stress;
                }
                else if (stress < min)
                {
                    min = stress;
                }
            }
            
            List<double> colorList = new List<double>();

            for (int i = 0; i < faces.Count; i++)
            {
                double stress = stresses[i*6+direction];

                R.Add(0);
                G.Add(0);
                B.Add(0);

                if (stress >= max*0.5 && max != 0)
                {
                    R[i] = 255;
                    G[i] = Convert.ToInt32(Math.Round(255 * (stress - max * 0.5) / (max * 0.5)));
                }
                else if (stress < max*0.5 && stress >= 0 && max != 0)
                {
                    G[i] = 255;
                    R[i] = Convert.ToInt32(Math.Round(255 * (stress) / (max * 0.5)));
                }
                else if (stress < 0 && stress > min*0.5 && min != 0)
                {
                    G[i] = 255;
                    B[i] = Convert.ToInt32(Math.Round(255 * (stress) / (min * 0.5)));
                }
                else if (stress <= min*0.5 && min != 0)
                {
                    B[i] = 255;
                    G[i] = Convert.ToInt32(Math.Round(255 * (stress - min*0.5) / (min * 0.5)));
                }
            }

            for (int i = 0; i < vertices.Count; i++)
            {
                List<int> vertex = new List<int>();
                int vR = 0, vG = 0, vB = 0;
                for (int j = 0; j < faces.Count; j++)
                {
                    if (faces[j].A == i || faces[j].B == i || faces[j].C == i)
                    {
                        vertex.Add(j);
                    }
                }
                for (int j = 0; j < vertex.Count; j++)
                {
                    vR += R[vertex[j]];
                    vG += G[vertex[j]];
                    vB += B[vertex[j]];
                }
                vR /= vertex.Count;
                vG /= vertex.Count;
                vB /= vertex.Count;

                meshOut.VertexColors.Add(vR, vG, vB);
            }
        }

        private List<Point3d> CreatePointList(List<Line> geometry)
        {
            List<Point3d> points = new List<Point3d>();

            for (int i = 0; i < geometry.Count; i++) //adds every point unless it already exists in list
            {
                Line l1 = geometry[i];
                if (!points.Contains(l1.From))
                {
                    points.Add(l1.From);
                }
                if (!points.Contains(l1.To))
                {
                    points.Add(l1.To);
                }
            }

            return points;
        }

        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                return Properties.Resources.Draw;
            }
        }

        public override Guid ComponentGuid
        {
            get { return new Guid("4b28fb40-2e66-4d19-a629-c630c079725a"); }
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

                Rectangle rec3 = rec2;
                rec3.X = rec2.Right + 2;

                Bounds = rec0;
                ButtonBounds = rec1;

            }

            GH_Palette xColor = GH_Palette.Black;

            private Rectangle ButtonBounds { get; set; }

            protected override void Render(GH_Canvas canvas, Graphics graphics, GH_CanvasChannel channel)
            {
                base.Render(canvas, graphics, channel);
                if (channel == GH_CanvasChannel.Objects)
                {
                    GH_Capsule button;
                    if (startDef == false)
                    {
                        button = GH_Capsule.CreateTextCapsule(ButtonBounds, ButtonBounds, GH_Palette.Grey, "Hidden", 3, 0);
                    }
                    else
                    {
                        button = GH_Capsule.CreateTextCapsule(ButtonBounds, ButtonBounds, xColor, "Displayed", 3, 0);
                    }
                    button.Render(graphics, Selected, false, false);
                    button.Dispose();
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
                        if (xColor == GH_Palette.Black) { DeformedGeometry.setToggles("Run", true); }
                        if (xColor == GH_Palette.Grey) { DeformedGeometry.setToggles("Run", false); }
                        sender.Refresh();
                        Owner.ExpireSolution(true);
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
            }
        }
    }



}