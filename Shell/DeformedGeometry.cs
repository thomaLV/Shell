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
          : base("DeformedShell", "DefS",
              "Displays the deformed shell, with or without coloring",
              "Koala", "Shell")
        {
        }

        //Initialize startcondition and polynomial order
        static bool startDef = true;
        static bool setColor = false;
        static bool X = false;
        static bool Y = false;
        static bool VonMisesButton = false;
        static bool RX = false;
        static bool RY = false;

        //Method to allow c hanging of variables via GUI (see Component Visual)
        public static void setToggles(string s, bool i)
        {
            if (s == "Run")
            {
                startDef = i;
            }
            if (s == "setColor")
            {
                setColor = i;
            }
            if (s == "X")
            {
                X = i;
            }
            if (s == "Y")
            {
                Y = i;
            }
            if (s == "VonMises")
            {
                VonMisesButton = i;
            }
            if (s == "RX")
            {
                RX= i;
            }
            if (s == "RY")
            {
                RY = i;
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
            pManager.AddMeshParameter("Mesh", "M", "Input Geometry (Mesh format)", GH_ParamAccess.item);
            pManager.AddNumberParameter("Scale", "S", "The Scale Factor for Deformation", GH_ParamAccess.item, 10);
            pManager.AddNumberParameter("Yield Strength", "YieldS", "The Yield Strength in MPa", GH_ParamAccess.list, new List<double> { 0, 0 });
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
                List<double> yieldStrength = new List<double>();
                List<Line> defGeometry = new List<Line>();
                List<Point3d> defPoints = new List<Point3d>();

                int[] h = new int[] { 0, 0, 0 };
                //Set expected inputs from Indata
                if (!DA.GetDataList(0, def)) return;
                if (!DA.GetDataList(1, stresses)) return;
                if (!DA.GetData(2, ref mesh)) return;
                if (!DA.GetData(3, ref scale)) return;
                if (!DA.GetDataList(4, yieldStrength)) return;
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

                int dimension = 123;
                if (X)
                {
                    dimension = 0;
                }
                else if (Y)
                {
                    dimension = 1;
                }
                else if (VonMisesButton)
                {
                    dimension = 7;
                    #region Von Mises
                    for (int j = 0; j < faces.Count; j++)
                    {
                        double sigma11 = stresses[j*6];
                        if (sigma11 >= 0)
                        {
                            sigma11 += Math.Abs(stresses[j * 6 + 3]);
                        }
                        else
                        {
                            sigma11 += -Math.Abs(stresses[j * 6 + 3]);
                        }

                        double sigma22 = stresses[j*6+1];
                        if (sigma22 >= 0)
                        {
                            sigma22 += Math.Abs(stresses[j * 6 + 4]);
                        }
                        else
                        {
                            sigma22 += -Math.Abs(stresses[j * 6 + 4]);
                        }

                        double sigma12 = stresses[j*6+2];
                        if (sigma12 >= 0)
                        {
                            sigma12 += Math.Abs(stresses[j * 6 + 5]);
                        }
                        else
                        {
                            sigma12 += -Math.Abs(stresses[j * 6 + 5]);
                        }

                        VonMises.Add(Math.Sqrt(sigma11 * sigma11 - sigma11 * sigma22 + sigma22 * sigma22 + 3 * sigma12 * sigma12));
                    }
                    
                    #endregion
                }
                else if (RX)
                {
                    dimension = 3;
                }
                else if (RY)
                {
                    dimension = 4;
                }

                Mesh coloredDefMesh = defmesh.DuplicateMesh();
                if (setColor && (stresses.Count > 1 || (stresses.Count == 1 && stresses[0] != 0) || VonMises.Count > 1 || (VonMises.Count == 1 && VonMises[0] != 0)) && (dimension < 8))
                {
                    // Direction can be 0:x
                    SetMeshColors(defmesh, stresses, VonMises, new_vertices, faces, dimension, yieldStrength, out coloredDefMesh);
                }
                 
                
                //Set output data
                DA.SetData(0, coloredDefMesh);
            }
        }   //End of main program

        private void SetMeshColors(Mesh meshIn, List<double> stresses, List<double> VonMises, List<Point3d> vertices, List<MeshFace> faces, int direction, List<double> yieldStrength, out Mesh meshOut)
        {
            meshOut = meshIn.DuplicateMesh();

            List<int> R = new List<int>(faces.Count);
            List<int> G = new List<int>(faces.Count);
            List<int> B = new List<int>(faces.Count);
            int[,] facesConnectedToVertex = new int[faces.Count,3];

            double max = 0;
            double min = 0;

            if (yieldStrength.Count == 1 && yieldStrength[0] > 1)
            {
                max = yieldStrength[0];
                min = -yieldStrength[0];
            }
            else if ((yieldStrength.Count == 1 && yieldStrength[0] == 0) || (yieldStrength[0] == 0 && yieldStrength[1] == 0) || yieldStrength.Count == 0)
            {
                for (int i = 0; i < stresses.Count / 6; i++)
                    {
                        double stress;
                        if (direction < 6)
                        {
                            stress = stresses[i * 6 + direction];
                        }
                        else
                        {
                            stress = VonMises[i];
                        }
                        if (stress > max)
                        {
                            max = stress;
                    }
                    else if (stress < min)
                    {
                        min = stress;
                    }
                }
            }
            else
            {
                if (yieldStrength[0] >= 0 && yieldStrength[1] <= 0)
                {
                    max = yieldStrength[0];
                    min = yieldStrength[1];
                }
                else if (yieldStrength[1] >= 0 && yieldStrength[0] <= 0)
                {
                    max = yieldStrength[1];
                    min = yieldStrength[0];
                }
                else
                {
                    AddRuntimeMessage(GH_RuntimeMessageLevel.Warning, "Warning message here");
                }
                
            }

            

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

                if (setColor)
                {
                    rec0.Height += 82;
                }
                else
                {
                    rec0.Height += 22;
                }

                Rectangle rec1 = rec0;
                rec1.X = rec0.Left + 1;
                
                if (setColor)
                {
                    rec1.Y = rec0.Bottom - 82;
                }
                else
                {
                    rec1.Y = rec0.Bottom - 22;
                }
                rec1.Width = (rec0.Width) / 2;
                rec1.Height = 22;
                rec1.Inflate(-2, -2);

                Rectangle rec2 = rec1;
                rec2.X = rec1.Right + 2;

                Rectangle rec3 = rec2;
                rec3.X = rec1.X;
                rec3.Y = rec1.Bottom + 2;

                Rectangle rec4 = rec3;
                rec4.X = rec3.X;
                rec4.Y = rec3.Bottom + 2;

                Rectangle rec5 = rec3;
                rec5.X = rec4.X;
                rec5.Y = rec4.Bottom + 2;

                Rectangle rec6 = rec3;
                rec6.X = rec5.Right + 2;
                rec6.Y = rec3.Bottom + 2;

                Rectangle rec7 = rec3;
                rec7.X = rec6.X;
                rec7.Y = rec6.Bottom + 2;

                Bounds = rec0;
                ButtonBounds = rec1;
                ButtonBounds1 = rec2;
                ButtonBounds2 = rec3;
                ButtonBounds3 = rec4;
                ButtonBounds4 = rec5;
                ButtonBounds5 = rec6;
                ButtonBounds6 = rec7;

            }

            GH_Palette displayed = GH_Palette.Black;
            GH_Palette setcolor = GH_Palette.Grey;
            GH_Palette xColor = GH_Palette.Grey;
            GH_Palette yColor = GH_Palette.Grey;
            GH_Palette VonMisesColor = GH_Palette.Grey;
            GH_Palette rxColor = GH_Palette.Grey;
            GH_Palette ryColor = GH_Palette.Grey;

            private Rectangle ButtonBounds { get; set; }
            private Rectangle ButtonBounds1 { get; set; }
            private Rectangle ButtonBounds2 { get; set; }
            private Rectangle ButtonBounds3 { get; set; }
            private Rectangle ButtonBounds4 { get; set; }
            private Rectangle ButtonBounds5 { get; set; }
            private Rectangle ButtonBounds6 { get; set; }

            protected override void Render(GH_Canvas canvas, Graphics graphics, GH_CanvasChannel channel)
            {
                base.Render(canvas, graphics, channel);
                if (channel == GH_CanvasChannel.Objects)
                {
                    GH_Capsule button;
                    if (startDef == false)
                    {
                        button = GH_Capsule.CreateTextCapsule(ButtonBounds, ButtonBounds, displayed, "Hidden", 3, 0);
                        button.Render(graphics, Selected, Owner.Locked, false);
                        button.Dispose();
                    }
                    else
                    {
                        button = GH_Capsule.CreateTextCapsule(ButtonBounds, ButtonBounds, displayed, "Displayed", 3, 0);
                        button.Render(graphics, Selected, Owner.Locked, false);
                        button.Dispose();
                    }
                    if (setColor == true)
                    {
                        GH_Capsule button2 = GH_Capsule.CreateTextCapsule(ButtonBounds1, ButtonBounds1, setcolor, "Colored", 2, 0);
                        button2.Render(graphics, Selected, Owner.Locked, false);
                        button2.Dispose();
                    }
                    else
                    {
                        GH_Capsule button2 = GH_Capsule.CreateTextCapsule(ButtonBounds1, ButtonBounds1, setcolor, "Uncolored", 2, 0);
                        button2.Render(graphics, Selected, Owner.Locked, false);
                        button2.Dispose();
                    }
                    if (setColor == true)
                    {
                        GH_Capsule button3 = GH_Capsule.CreateTextCapsule(ButtonBounds2, ButtonBounds2, xColor, "X Stresses", 2, 0);
                        button3.Render(graphics, Selected, Owner.Locked, false);
                        button3.Dispose();
                    }
                    if (setColor == true)
                    {
                        GH_Capsule button4 = GH_Capsule.CreateTextCapsule(ButtonBounds3, ButtonBounds3, yColor, "Y Stresses", 2, 0);
                        button4.Render(graphics, Selected, Owner.Locked, false);
                        button4.Dispose();
                    }
                    if (setColor == true)
                    {
                        GH_Capsule button5 = GH_Capsule.CreateTextCapsule(ButtonBounds4, ButtonBounds4, VonMisesColor, "Von Mises", 2, 0);
                        button5.Render(graphics, Selected, Owner.Locked, false);
                        button5.Dispose();
                    }
                    if (setColor == true)
                    {
                        GH_Capsule button6 = GH_Capsule.CreateTextCapsule(ButtonBounds5, ButtonBounds5, rxColor, "RX Stresses", 2, 0);
                        button6.Render(graphics, Selected, Owner.Locked, false);
                        button6.Dispose();
                    }
                    if (setColor == true)
                    {
                        GH_Capsule button7 = GH_Capsule.CreateTextCapsule(ButtonBounds6, ButtonBounds6, ryColor, "RY Stresses", 2, 0);
                        button7.Render(graphics, Selected, Owner.Locked, false);
                        button7.Dispose();
                    }
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
                    }
                    rec = ButtonBounds1;
                    if (rec.Contains(e.CanvasLocation))
                    {
                        switchColor("setColor");
                    }
                    rec = ButtonBounds2;
                    if (rec.Contains(e.CanvasLocation))
                    {
                        switchColor("X");
                    }
                    rec = ButtonBounds3;
                    if (rec.Contains(e.CanvasLocation))
                    {
                        switchColor("Y");
                    }
                    rec = ButtonBounds4;
                    if (rec.Contains(e.CanvasLocation))
                    {
                        switchColor("VonMises");
                    }
                    rec = ButtonBounds5;
                    if (rec.Contains(e.CanvasLocation))
                    {
                        switchColor("RX");
                    }
                    rec = ButtonBounds6;
                    if (rec.Contains(e.CanvasLocation))
                    {
                        switchColor("RY");
                    }

                    if (displayed == GH_Palette.Black) { DeformedGeometry.setToggles("Run", true); }
                    if (displayed == GH_Palette.Grey) { DeformedGeometry.setToggles("Run", false); }
                    if (setcolor == GH_Palette.Black) { DeformedGeometry.setToggles("setColor", true); }
                    if (setcolor == GH_Palette.Grey) { DeformedGeometry.setToggles("setColor", false); }
                    if (xColor == GH_Palette.Black) { DeformedGeometry.setToggles("X", true); }
                    if (xColor == GH_Palette.Grey) { DeformedGeometry.setToggles("X", false); }
                    if (yColor == GH_Palette.Black) { DeformedGeometry.setToggles("Y", true); }
                    if (yColor == GH_Palette.Grey) { DeformedGeometry.setToggles("Y", false); }
                    if (VonMisesColor == GH_Palette.Black) { DeformedGeometry.setToggles("VonMises", true); }
                    if (VonMisesColor == GH_Palette.Grey) { DeformedGeometry.setToggles("VonMises", false); }
                    if (rxColor == GH_Palette.Black) { DeformedGeometry.setToggles("RX", true); }
                    if (rxColor == GH_Palette.Grey) { DeformedGeometry.setToggles("RX", false); }
                    if (ryColor == GH_Palette.Black) { DeformedGeometry.setToggles("RY", true); }
                    if (ryColor == GH_Palette.Grey) { DeformedGeometry.setToggles("RY", false); }
                    sender.Refresh();
                    Owner.ExpireSolution(true);
                    return GH_ObjectResponse.Handled;
                    
                }
                return base.RespondToMouseDown(sender, e);
            }

            private void switchColor(string button)
            {
                if (button == "Run")
                {
                    if (displayed == GH_Palette.Black) { displayed = GH_Palette.Grey; }
                    else { displayed = GH_Palette.Black; }
                }
                if (button == "setColor")
                {
                    if (setcolor == GH_Palette.Black)
                    {
                        setcolor = GH_Palette.Grey;
                        xColor = GH_Palette.Grey;
                        yColor = GH_Palette.Grey;
                        VonMisesColor = GH_Palette.Grey;
                        rxColor = GH_Palette.Grey;
                        ryColor = GH_Palette.Grey;
                    }
                    else { setcolor = GH_Palette.Black; }
                }
                if (button == "X" && setcolor == GH_Palette.Black)
                {
                    if (xColor == GH_Palette.Black) { xColor = GH_Palette.Grey; }
                    else
                    {
                        xColor = GH_Palette.Black;
                        yColor = GH_Palette.Grey;
                        VonMisesColor = GH_Palette.Grey;
                        rxColor = GH_Palette.Grey;
                        ryColor = GH_Palette.Grey;
                    }
                }
                if (button == "Y" && setcolor == GH_Palette.Black)
                {
                    if (yColor == GH_Palette.Black) { yColor = GH_Palette.Grey; }
                    else
                    {
                        yColor = GH_Palette.Black;
                        xColor = GH_Palette.Grey;
                        VonMisesColor = GH_Palette.Grey;
                        rxColor = GH_Palette.Grey;
                        ryColor = GH_Palette.Grey;
                    }
                }
                if (button == "VonMises" && setcolor == GH_Palette.Black)
                {
                    if (VonMisesColor == GH_Palette.Black) { VonMisesColor = GH_Palette.Grey; }
                    else
                    {
                        VonMisesColor = GH_Palette.Black;
                        xColor = GH_Palette.Grey;
                        yColor = GH_Palette.Grey;
                        rxColor = GH_Palette.Grey;
                        ryColor = GH_Palette.Grey;
                    }
                }
                if (button == "RX" && setcolor == GH_Palette.Black)
                {
                    if (rxColor == GH_Palette.Black) { rxColor = GH_Palette.Grey; }
                    else
                    {
                        xColor = GH_Palette.Grey;
                        yColor = GH_Palette.Grey;
                        VonMisesColor = GH_Palette.Grey;
                        rxColor = GH_Palette.Black;
                        ryColor = GH_Palette.Grey;
                    }
                }
                if (button == "RY" && setcolor == GH_Palette.Black)
                {
                    if (ryColor == GH_Palette.Black) { ryColor = GH_Palette.Grey; }
                    else
                    {
                        xColor = GH_Palette.Grey;
                        yColor = GH_Palette.Grey;
                        VonMisesColor = GH_Palette.Grey;
                        rxColor = GH_Palette.Grey;
                        ryColor = GH_Palette.Black;
                    }
                }
            }
        }
    }



}