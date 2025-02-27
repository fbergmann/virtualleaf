/*
 *  This file is part of the VirtualLeaf.
 *
 *  VirtualLeaf is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  VirtualLeaf is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with the Virtual Leaf.  If not, see <http://www.gnu.org/licenses/>.
 *
 *  Copyright 2011-14 Dirk De Vos, Jan Broeckhove, UA/CoMP.
 */
/**
 * @file
 * Interface for Wortel plug-in.
 */
#include <string>
#include <QObject>
#include <QtGui>

#include "simplugin.h"

static const std::string Wortel_header_id("$Id$");

class Mesh;

class Wortel : public QObject, SimPluginInterface {
  Q_OBJECT
  Q_INTERFACES(SimPluginInterface)

public:
  Wortel();

  // see SimPluginInterface
  virtual QString ModelID();

  // default LeafML file
  virtual QString DefaultLeafML();

  // return number of chemicals
  virtual int NChem();

  /// Set mesh pointer at startup.
  virtual void SetMesh(Mesh* m);

  // Executed after the cellular mechanics steps have equillibrized
  virtual void CellHouseKeeping (CellBase* c);

  // Differential equations describing transport of chemicals from cell to cell
  virtual void CelltoCellTransport(Wall* w, double* dchem_c1, double* dchem_c2);

  // Differential equations describing chemical reactions taking place at or near the cell walls
  // (e.g. PIN accumulation)
  virtual void WallDynamics(Wall* w, double* dw1, double* dw2);

  // Differential equations describing chemical reactions inside the cells
  virtual void CellDynamics(CellBase* c, double* dchem);

  // to be executed after a cell division
  virtual void OnDivide(ParentInfo* parent_info, CellBase* daughter1, CellBase* daughter2);

  // to be executed for coloring a cell
  virtual void SetCellColor(CellBase* c, QColor* color);

  // see SimPluginInterface
  virtual double PINflux(CellBase* this_cell, CellBase* adjacent_cell, Wall* w);

private:
  /// Points to the mesh this model works with.
  Mesh*   m_mesh;
  double time_now;
  double time_start;
  int mModel;

  // Calculate rate according to Michaelis Menten kinetics
  double Michaelis(double M1, double J1, double K1, double S1);

  // Calculate rate according to Goldbeter-Koshland kinetics
  double Goldbeter(double A1, double A2, double A3, double A4);

  //Calculate rate according to Hill cooperative kinetics
  double Hill(double Vm, double Km, double h, double S);

};




