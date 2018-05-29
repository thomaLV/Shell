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
            pManager.AddNumberParameter("Stresses", "Stress", "Stresses from ShellCalc", GH_ParamAccess.list, new List<double> { 0 });
            pManager.AddNumberParameter("VonMises Stress", "VM", "VonMises from ShellCalc", GH_ParamAccess.list, new List<double> { 0 });
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
                List<double> VonMises = new List<double>();
                Mesh mesh = new Mesh();
                double scale = 10;
                List<Line> defGeometry = new List<Line>();
                List<Point3d> defPoints = new List<Point3d>();

                int dimension = 0;
                int[] h = new int[] { 0, 0, 0 };
                //Set expected inputs from Indata
                if (!DA.GetDataList(0, def)) return;
                if (!DA.GetDataList(1, stresses)) return;
                if (!DA.GetDataList(2, VonMises)) return;
                if (!DA.GetData(3, ref mesh)) return;
                if (!DA.GetData(4, ref scale)) return;
                if (!DA.GetData(5, ref dimension)) return;
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
                #endregion

                Mesh coloredDefMesh = defmesh.DuplicateMesh();
                if (stresses.Count > 1 || (stresses.Count == 1 && stresses[0] != 0) || VonMises.Count > 1 || (VonMises.Count == 1 && VonMises[0] != 0))
                {
                    // Direction can be 0:x
                    SetMeshColors(defmesh, stresses, VonMises, new_vertices, faces, dimension, out coloredDefMesh);
                }
                 
                
                //Set output data
                DA.SetData(0, coloredDefMesh);
            }
        }   //End of main program

        private void SetMeshColors(Mesh meshIn, List<double> stresses, List<double> VonMises, List<Point3d> vertices, List<MeshFace> faces, int direction, out Mesh meshOut)
        {
            meshOut = meshIn.DuplicateMesh();

            double max = 355;
            double min = -355;
            List<int> R = new List<int>(faces.Count);
            List<int> G = new List<int>(faces.Count);
            List<int> B = new List<int>(faces.Count);
            int[,] facesConnectedToVertex = new int[faces.Count,3];

            //for (int i = 0; i < stresses.Count / 6; i++)
            //{
            //    double stress = stresses[i * 6 + direction];
            //    if (stress > max)
            //    {
            //        max = stress;
            //    }
            //    else if (stress < min)
            //    {
            //        min = stress;
            //    }
            //}

            List<double> colorList = new List<double>();

            for (int i = 0; i < faces.Count; i++)
            {
                double stress;
                if (direction < 6)
                {
                    stress = stresses[i*6+direction];
                }
                else
                {
                    stress = VonMises[i];
                }
                

                R.Add(0);
                G.Add(0);
                B.Add(0);

                if (stress >= max)
                {
                    R[i] = 255;
                }
                else if (stress >= max*0.5 && max != 0)
                {
                    R[i] = 255;
                    G[i] = Convert.ToInt32(Math.Round(255 * (1 - (stress - max * 0.5) / (max * 0.5))));
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
                else if (stress <= min*0.5 && min != 0 && stress > min)
                {
                    B[i] = 255;
                    G[i] = Convert.ToInt32(Math.Round(255 * (1 - (stress - min * 0.5) / (min * 0.5))));
                }
                else if (stress <= min)
                {
                    B[i] = 255;
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