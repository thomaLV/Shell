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
    public class BDCComponent : GH_Component
    {
        public BDCComponent()
          : base("BDCComponent", "BDCs",
              "Description",
              "Koala", "Shell")
        {
        }

        //Initialize BDCs
        static int x = 0;
        static int y = 0;
        static int z = 0;
        static int rx = 0;
        static int ry = 0;
        static int rz = 0;


        //Method to allow c hanging of variables via GUI (see Component Visual)
        public static void setBDC(string s, int i)
        {
            if (s == "X")
            {
                x = i;
            }
            else if (s == "Y")
            {
                y = i;
            }
            else if (s == "Z")
            {
                z = i;
            }
            else if (s == "RX")
            {
                rx = i;
            }
        }

        public override void CreateAttributes()
        {
            m_attributes = new Attributes_Custom(this);
        }

        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddPointParameter("Points", "P", "Points to apply Boundary Conditions", GH_ParamAccess.list);
        }

        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddTextParameter("B.Cond.", "BDC", "Boundary Conditions for 3D Beam Calculation", GH_ParamAccess.list);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            #region Fetch inputs
            //Expected inputs
            List<Point3d> pointList = new List<Point3d>();          //List of points where BDC is to be applied
            List<string> pointInStringFormat = new List<string>();  //output in form of list of strings


            //Set expected inputs from Indata and aborts with error message if input is incorrect
            if (!DA.GetDataList(0, pointList)) return;
            #endregion

            #region Format output
            string BDCString = x + "," + y + "," + z + "," + rx + "," + ry + "," + rz;

            for (int i = 0; i < pointList.Count; i++)   //Format stringline for all points (identical boundary conditions for all points)
            {
                pointInStringFormat.Add(pointList[i].X + "," + pointList[i].Y + "," + pointList[i].Z + ":" + BDCString);
            }
            #endregion

            DA.SetDataList(0, pointInStringFormat);
        } //End of main program

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
                return Properties.Resources.BDCs;
            }
        }

        public override Guid ComponentGuid
        {
            get { return new Guid("58ccdcb8-b1c3-411b-b501-c91a46665e86"); }
        }

        /// Component Visual//
        public class Attributes_Custom : Grasshopper.Kernel.Attributes.GH_ComponentAttributes
        {
            public Attributes_Custom(GH_Component owner) : base(owner) { }
            protected override void Layout()
            {
                base.Layout();

                Rectangle rec0 = GH_Convert.ToRectangle(Bounds);

                rec0.Height += 42;

                Rectangle rec1 = rec0;
                rec1.X = rec0.Left + 1;
                rec1.Y = rec0.Bottom - 42;
                rec1.Width = (rec0.Width) / 3 + 1;
                rec1.Height = 22;
                rec1.Inflate(-2, -2);

                Rectangle rec2 = rec1;
                rec2.X = rec1.Right + 2;

                Rectangle rec3 = rec2;
                rec3.X = rec2.Right + 2;

                Rectangle rec4 = rec1;
                rec4.Y = rec1.Bottom + 2;

                Rectangle rec5 = rec4;
                rec5.X = rec4.Right + 2;

                Rectangle rec6 = rec5;
                rec6.X = rec2.Right + 2;

                Bounds = rec0;
                BoundsAllButtons = rec0;
                ButtonBounds = rec1;
                ButtonBounds2 = rec2;
                ButtonBounds3 = rec3;
                ButtonBounds4 = rec4;
                ButtonBounds5 = rec5;
                ButtonBounds6 = rec6;

            }

            GH_Palette xColor = GH_Palette.Black;
            GH_Palette yColor = GH_Palette.Black;
            GH_Palette zColor = GH_Palette.Black;
            GH_Palette rxColor = GH_Palette.Black;
            GH_Palette ryColor = GH_Palette.Black;
            GH_Palette rzColor = GH_Palette.Black;
            
            private Rectangle BoundsAllButtons { get; set; }
            private Rectangle ButtonBounds { get; set; }
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
                    GH_Capsule button = GH_Capsule.CreateTextCapsule(ButtonBounds, ButtonBounds, xColor, "X", 3, 0);
                    button.Render(graphics, Selected, false, false);
                    button.Dispose();
                }
                if (channel == GH_CanvasChannel.Objects)
                {
                    GH_Capsule button2 = GH_Capsule.CreateTextCapsule(ButtonBounds2, ButtonBounds2, yColor, "Y", 2, 0);
                    button2.Render(graphics, Selected, Owner.Locked, false);
                    button2.Dispose();
                }
                if (channel == GH_CanvasChannel.Objects)
                {
                    GH_Capsule button3 = GH_Capsule.CreateTextCapsule(ButtonBounds3, ButtonBounds3, zColor, "Z", 2, 0);
                    button3.Render(graphics, Selected, Owner.Locked, false);
                    button3.Dispose();
                }
                if (channel == GH_CanvasChannel.Objects)
                {
                    GH_Capsule button4 = GH_Capsule.CreateTextCapsule(ButtonBounds4, ButtonBounds4, rxColor, "RX", 2, 0);
                    button4.Render(graphics, Selected, Owner.Locked, false);
                    button4.Dispose();
                }
                if (channel == GH_CanvasChannel.Objects)
                {
                    GH_Capsule button5 = GH_Capsule.CreateTextCapsule(ButtonBounds5, ButtonBounds5, ryColor, "RY", 2, 0);
                    button5.Render(graphics, Selected, Owner.Locked, false);
                    button5.Dispose();
                }
                if (channel == GH_CanvasChannel.Objects)
                {
                    GH_Capsule button6 = GH_Capsule.CreateTextCapsule(ButtonBounds6, ButtonBounds6, rzColor, "RZ", 2, 0);
                    button6.Render(graphics, Selected, Owner.Locked, false);
                    button6.Dispose();
                }
            }

            public override GH_ObjectResponse RespondToMouseDown(GH_Canvas sender, GH_CanvasMouseEvent e)
            {
                if (e.Button == MouseButtons.Left)
                {
                    RectangleF rec = ButtonBounds;
                    if (rec.Contains(e.CanvasLocation))
                    {
                        switchColor("X");
                    }
                    rec = ButtonBounds2;
                    if (rec.Contains(e.CanvasLocation))
                    {
                        switchColor("Y");
                    }
                    rec = ButtonBounds3;
                    if (rec.Contains(e.CanvasLocation))
                    {
                        switchColor("Z");
                    }
                    rec = ButtonBounds4;
                    if (rec.Contains(e.CanvasLocation))
                    {
                        switchColor("RX");
                    }
                    rec = ButtonBounds5;
                    if (rec.Contains(e.CanvasLocation))
                    {
                        switchColor("RY");
                    }
                    rec = ButtonBounds6;
                    if (rec.Contains(e.CanvasLocation))
                    {
                        switchColor("RZ");
                    }
                    rec = BoundsAllButtons;
                    if (rec.Contains(e.CanvasLocation))
                    {
                        if (xColor == GH_Palette.Black) { BDCComponent.setBDC("X", 0); }
                        if (xColor == GH_Palette.Grey) { BDCComponent.setBDC("X", 1); }
                        if (yColor == GH_Palette.Black) { BDCComponent.setBDC("Y", 0); }
                        if (yColor == GH_Palette.Grey) { BDCComponent.setBDC("Y", 1); }
                        if (zColor == GH_Palette.Black) { BDCComponent.setBDC("Z", 0); }
                        if (zColor == GH_Palette.Grey) { BDCComponent.setBDC("Z", 1); }
                        if (rxColor == GH_Palette.Black) { BDCComponent.setBDC("RX", 0); }
                        if (rxColor == GH_Palette.Grey) { BDCComponent.setBDC("RX", 1); }
                        if (ryColor == GH_Palette.Black) { BDCComponent.setBDC("RY", 0); }
                        if (ryColor == GH_Palette.Grey) { BDCComponent.setBDC("RY", 1); }
                        if (rzColor == GH_Palette.Black) { BDCComponent.setBDC("RZ", 0); }
                        if (rzColor == GH_Palette.Grey) { BDCComponent.setBDC("RZ", 1); }
                        sender.Refresh();
                        Owner.ExpireSolution(true);
                    }
                    return GH_ObjectResponse.Handled;
                }
                return base.RespondToMouseDown(sender, e);
            }

            private void switchColor(string button)
            {
                if (button == "X")
                {
                    if (xColor == GH_Palette.Black) { xColor = GH_Palette.Grey; }
                    else { xColor = GH_Palette.Black; }
                }
                else if (button == "Y")
                {
                    if (yColor == GH_Palette.Black) { yColor = GH_Palette.Grey; }
                    else { yColor = GH_Palette.Black; }
                }
                else if (button == "Z")
                {
                    if (zColor == GH_Palette.Black) { zColor = GH_Palette.Grey; }
                    else { zColor = GH_Palette.Black; }
                }
                else if (button == "RX")
                {
                    if (rxColor == GH_Palette.Black) { rxColor = GH_Palette.Grey; }
                    else { rxColor = GH_Palette.Black; }
                }
                else if (button == "RY")
                {
                    if (ryColor == GH_Palette.Black) { ryColor = GH_Palette.Grey; }
                    else { ryColor = GH_Palette.Black; }
                }
                else if (button == "RZ")
                {
                    if (rzColor == GH_Palette.Black) { rzColor = GH_Palette.Grey; }
                    else { rzColor = GH_Palette.Black; }
                }
            }
        }
    }
}
