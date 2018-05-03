using System;
using System.Collections.Generic;

using Grasshopper.Kernel;
using Rhino.Geometry;
using System.Drawing;
using Grasshopper.GUI.Canvas;
using System.Windows.Forms;
using Grasshopper.GUI;


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

        //Initialize startcondition
        static bool startDef = true;
        static bool p2 = true;
        static bool p3 = false;


        //Method to allow c hanging of variables via GUI (see Component Visual)
        public static void setToggles(string s, bool i)
        {
            if (s == "Run")
            {
                startDef = i;
            }
            else if (s == "2nd")
            {
                p2 = i;
            }
            else if (s == "3rd")
            {
                p3 = i;
            }
            Grasshopper.Instances.ActiveCanvas.Document.ExpireSolution();
            Grasshopper.Instances.ActiveCanvas.Document.NewSolution(false);
        }

        public override void CreateAttributes()
        {
            m_attributes = new Attributes_Custom(this);
        }

        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddNumberParameter("Deformation", "Def", "Deformations from 3DBeamCalc", GH_ParamAccess.list);
            pManager.AddLineParameter("Geometry", "G", "Input Geometry (Line format)", GH_ParamAccess.list);
            pManager.AddNumberParameter("Scale", "S", "The Scale Factor for Deformation", GH_ParamAccess.item, 1000);
        }

        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddLineParameter("Deformed Geometry", "Def.G.", "Deformed Geometry as List of Lines", GH_ParamAccess.list);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            if (startDef)
            {
                #region Fetch input
                //Expected inputs and outputs
                List<double> def = new List<double>();
                List<Line> geometry = new List<Line>();
                double scale = 1000;
                List<Line> defGeometry = new List<Line>();
                List<Point3d> defPoints = new List<Point3d>();

                //Set expected inputs from Indata
                if (!DA.GetDataList(0, def)) return;
                if (!DA.GetDataList(1, geometry)) return;
                if (!DA.GetData(2, ref scale)) return;
                #endregion

                #region Scale deformations
                //List all nodes (every node only once), numbering them according to list index
                List<Point3d> points = CreatePointList(geometry);

                int index = 0;
                //loops through all points and scales x-, y- and z-dir

                //u(x) = Na, N = shape func, a = nodal values (dof) 
                foreach (Point3d point in points)
                {
                    //fetch global x,y,z placement of point
                    double x = point.X;
                    double y = point.Y;
                    double z = point.Z;

                    //scales x and z according to input Scale
                    defPoints.Add(new Point3d(x + scale * def[index], y + scale * def[index + 1], z + scale * def[index + 2]));
                    index += 6;
                }
                #endregion

                #region Create geometry
                //creates deformed geometry based on initial geometry placement
                foreach (Line line in geometry)
                {
                    //fetches index of original start and endpoint
                    int i1 = points.IndexOf(line.From);
                    int i2 = points.IndexOf(line.To);

                    //creates new line based on scaled deformation of said points
                    defGeometry.Add(new Line(defPoints[i1], defPoints[i2]));
                }
                #endregion


                //Set output data
                DA.SetDataList(0, defGeometry);
            }
        }   //End of main program


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


        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
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
                ButtonBounds2 = rec2;
                ButtonBounds3 = rec3;

            }

            GH_Palette xColor = GH_Palette.Black;
            GH_Palette yColor = GH_Palette.Black;
            GH_Palette zColor = GH_Palette.Grey;

            private Rectangle ButtonBounds { get; set; }
            private Rectangle ButtonBounds2 { get; set; }
            private Rectangle ButtonBounds3 { get; set; }

            protected override void Render(GH_Canvas canvas, Graphics graphics, GH_CanvasChannel channel)
            {
                base.Render(canvas, graphics, channel);
                if (channel == GH_CanvasChannel.Objects)
                {
                    GH_Capsule button = GH_Capsule.CreateTextCapsule(ButtonBounds, ButtonBounds, xColor, "Run", 3, 0);
                    button.Render(graphics, Selected, false, false);
                    button.Dispose();
                }
                if (channel == GH_CanvasChannel.Objects)
                {
                    GH_Capsule button2 = GH_Capsule.CreateTextCapsule(ButtonBounds2, ButtonBounds2, yColor, "2nd", 2, 0);
                    button2.Render(graphics, Selected, Owner.Locked, false);
                    button2.Dispose();
                }
                if (channel == GH_CanvasChannel.Objects)
                {
                    GH_Capsule button3 = GH_Capsule.CreateTextCapsule(ButtonBounds3, ButtonBounds3, zColor, "3rd", 2, 0);
                    button3.Render(graphics, Selected, Owner.Locked, false);
                    button3.Dispose();
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
                        return GH_ObjectResponse.Handled;
                    }
                    rec = ButtonBounds2;
                    if (rec.Contains(e.CanvasLocation))
                    {
                        if (yColor == GH_Palette.Black) { DeformedGeometry.setToggles("2nd", false); DeformedGeometry.setToggles("3rd", true); }
                        if (yColor == GH_Palette.Grey) { DeformedGeometry.setToggles("2nd", true); DeformedGeometry.setToggles("3rd", false); }
                        switchColor("2nd");
                        switchColor("3rd");
                        sender.Refresh();
                        return GH_ObjectResponse.Handled;
                    }
                    rec = ButtonBounds3;
                    if (rec.Contains(e.CanvasLocation))
                    {
                        if (zColor == GH_Palette.Black) { DeformedGeometry.setToggles("3rd", false); DeformedGeometry.setToggles("2nd", false); }
                        if (zColor == GH_Palette.Grey) { DeformedGeometry.setToggles("3rd", true); DeformedGeometry.setToggles("2nd", true); }
                        switchColor("3rd");
                        switchColor("2nd");
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
                else if (button == "2nd")
                {
                    if (yColor == GH_Palette.Black) { yColor = GH_Palette.Grey; }
                    else { yColor = GH_Palette.Black; }
                }
                else if (button == "3rd")
                {
                    if (zColor == GH_Palette.Black) { zColor = GH_Palette.Grey; }
                    else { zColor = GH_Palette.Black; }
                }
            }
        }
    }



}