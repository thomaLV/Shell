using System;
using System.Collections.Generic;

using Grasshopper.Kernel;
using Rhino.Geometry;

namespace Shell
{
    public class SetLoads : GH_Component
    {
        public SetLoads()
          : base("PointLoads Shell", "PL",
              "Point loads to apply to a shell structure",
              "Koala", "Shell")
        {
        }
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddPointParameter("Points", "P", "Points to apply load(s)", GH_ParamAccess.list);
            pManager.AddNumberParameter("Load", "L", "Load originally given i Newtons (N), give one load for all points or list of loads for each point", GH_ParamAccess.list);
            pManager.AddNumberParameter("angle (xz)", "axz", "give angle for load in xz plane", GH_ParamAccess.list, 90);
            pManager.AddNumberParameter("angle (xy)", "axy", "give angle for load in xy plane", GH_ParamAccess.list, 0);
            //pManager[2].Optional = true; //Code can run without a given angle (90 degrees is initial value)
        }

        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddTextParameter("PointLoads", "PL", "PointLoads formatted for Calculation Component", GH_ParamAccess.list);
        }

        protected override void SolveInstance(IGH_DataAccess DA)
        {
            #region Fetch inputs
            //Expected inputs and output
            List<Point3d> pointList = new List<Point3d>();              //List of points where load will be applied
            List<double> loadList = new List<double>();                 //List or value of load applied
            List<double> anglexz = new List<double>();                  //Initial xz angle 90, angle from x axis in xz plane for load
            List<double> anglexy = new List<double>();                  //Initial xy angle 0, angle from x axis in xy plane for load
            List<string> pointInStringFormat = new List<string>();      //preallocate final string output

            //Set expected inputs from Indata
            if (!DA.GetDataList(0, pointList)) return;
            if (!DA.GetDataList(1, loadList)) return;
            if (!DA.GetDataList(2, anglexz)) return;
            if (!DA.GetDataList(3, anglexy)) return;
            #endregion

            #region Format pointloads
            //initialize temporary stringline and load vectors
            string vectorString;
            double load = 0;
            double xvec = 0;
            double yvec = 0;
            double zvec = 0;

            if (loadList.Count == 1 && anglexz.Count == 1)              //loads and angles are identical for all points 
            {
                load = -1 * loadList[0];                                //negativ load for z-dir
                xvec = Math.Round(load * Math.Cos(anglexz[0] * Math.PI / 180) * Math.Cos(anglexy[0] * Math.PI / 180), 5);
                yvec = Math.Round(load * Math.Cos(anglexz[0] * Math.PI / 180) * Math.Sin(anglexy[0] * Math.PI / 180), 5);
                zvec = Math.Round(load * Math.Sin(anglexz[0] * Math.PI / 180), 5);

                vectorString = xvec + "," + yvec + "," + zvec;
                for (int i = 0; i < pointList.Count; i++)               //adds identical load to all points in pointList
                {
                    pointInStringFormat.Add(pointList[i].X + "," + pointList[i].Y + "," + pointList[i].Z + ":" + vectorString);
                }
            }
            else   //loads and angles may be different => calculate new xvec, yvec, zvec for all loads                          
            {
                for (int i = 0; i < pointList.Count; i++)
                {
                    if (loadList.Count < i)             //if pointlist is larger than loadlist, set last load value in remaining points
                    {
                        vectorString = xvec + "," + yvec + "," + zvec;
                    }
                    else
                    {
                        load = -1 * loadList[i];        //negative load for z-dir

                        xvec = Math.Round(load * Math.Cos(anglexz[i]) * Math.Cos(anglexy[i]), 2);
                        yvec = Math.Round(load * Math.Cos(anglexz[i]) * Math.Sin(anglexy[i]), 2);
                        zvec = Math.Round(load * Math.Sin(anglexz[i]), 2);

                        vectorString = xvec + "," + yvec + "," + zvec;
                    }

                    pointInStringFormat.Add(pointList[i].X + "," + pointList[i].Y + "," + pointList[i].Z + ":" + vectorString);
                }
            }
            #endregion

            //Set output data
            DA.SetDataList(0, pointInStringFormat);
        }

        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                return Properties.Resources.Pointloads;
            }
        }

        public override Guid ComponentGuid
        {
            get { return new Guid("2935c931-2647-4bc5-b851-68e7d4af9001"); }
        }
    }
}