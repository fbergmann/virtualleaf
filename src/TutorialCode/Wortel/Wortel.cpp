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
 * Implementation for Wortel plug-in.
 */

#include <cmath>
#include "mesh.h"
#include "cellbase.h"
#include "wallbase.h"
#include "Wortel.h"
#include "parameter.h"

static const std::string Wortel_module_id("$Id$");
using std::string;
using std::pow;

Wortel::Wortel()
  : QObject()
  , m_mesh(NULL)
  , time_now(0)
  , time_start(0)
  , mModel(1) // valid values are 1 .. 12
{

}

QString Wortel::ModelID()
{
  return QString("Wortel");
}

QString Wortel::DefaultLeafML()
{
  if (mModel == 9)
    return QString("Wortel_2.xml");
  if (mModel == 12)
    return QString("Wortel_3.xml");
  return QString("Wortel_1.xml");
}
int Wortel::NChem()
{
  return 9;
}


void Wortel::SetMesh(Mesh* mesh)
{
  m_mesh = mesh;
}


void Wortel::OnDivide(ParentInfo* parent_info, CellBase* daughter1, CellBase* daughter2)
{

  // Species distributes between parent and daughter according to area
  double area1 = daughter1->Area(), area2 = daughter2->Area();
  double tot_area = area1 + area2;

  daughter1->SetChemical(0, daughter1->Chemical(0) * (area1 / tot_area));
  daughter2->SetChemical(0, daughter2->Chemical(0) * (area2 / tot_area));
  daughter1->SetChemical(1, daughter1->Chemical(1) * (area1 / tot_area));
  daughter2->SetChemical(1, daughter2->Chemical(1) * (area2 / tot_area));
  daughter1->SetChemical(2, daughter1->Chemical(2) * (area1 / tot_area));
  daughter2->SetChemical(2, daughter2->Chemical(2) * (area2 / tot_area));
  daughter1->SetChemical(3, daughter1->Chemical(3) * (area1 / tot_area));
  daughter2->SetChemical(3, daughter2->Chemical(3) * (area2 / tot_area));
  daughter1->SetChemical(4, daughter1->Chemical(4) * (area1 / tot_area));
  daughter2->SetChemical(4, daughter2->Chemical(4) * (area2 / tot_area));
  daughter1->SetChemical(5, daughter1->Chemical(5) * (area1 / tot_area));
  daughter2->SetChemical(5, daughter2->Chemical(5) * (area2 / tot_area));
  daughter1->SetChemical(6, daughter1->Chemical(6) * (area1 / tot_area));
  daughter2->SetChemical(6, daughter2->Chemical(6) * (area2 / tot_area));
  daughter1->SetChemical(7, daughter1->Chemical(7) * (area1 / tot_area));
  daughter2->SetChemical(7, daughter2->Chemical(7) * (area2 / tot_area));
  daughter1->SetChemical(8, daughter1->Chemical(8) * (area1 / tot_area));
  daughter2->SetChemical(8, daughter2->Chemical(8) * (area2 / tot_area));


  //	daughter1->SetTransporters(1);
  //	daughter2->SetTransporters(1);

  daughter1->SetTransporters(1, 1., 1);
  daughter2->SetTransporters(1, 1., 1);

  daughter1->SetAStrain(0.);
  daughter2->SetAStrain(0.);

  daughter1->SetPrevArea(fabs(area1));
  daughter2->SetPrevArea(fabs(area2));

  daughter1->SetCellCycleTime( m_mesh->getTime() ); //Records simulation time upon division
  daughter2->SetCellCycleTime( m_mesh->getTime() ); //id.

}


void Wortel::SetCellColor(CellBase* c, QColor* color)
{
  if ((c->Chemical(0)/( c->Area() )) <= 0.)
  {
    color->setRgb(255.0,0,0);
  }
  else
  {
    color->setRgb(255.0, ((c->Chemical(0)/( c->Area() )) / (1000 + ( c->Chemical(0) / ( c->Area() ) ))) * 255., 0);
  }
}


void Wortel::CellHouseKeeping(CellBase* c)
{

  //Calculate areal strain rate over rdt
  double Acurr = fabs( c->Area() );
  double Aprev = fabs( c->GetPrevArea());
  c->SetAStrain(( Acurr - Aprev ) / Aprev );
  c->SetPrevArea( Acurr);


  //Keep track of tip position
  static double tip_position = 660.;
  static double next_tip_position = 660.;
  if (c->Index() == 0)
  {
    tip_position = next_tip_position;
  }
  if (c->Index() == 194)
  {
    next_tip_position = c->Centroid().y;
  }


  ///////////////////// **GROWTH RULES** /////////////////////////


  // MODEL 12 (CF. TABLE S1 FOR PARAMETERS)

  if ( mModel == 12 &&  c->CellType() != 0 && ( c->Chemical(3) / c->Area() ) >= 0.7 )
  {
    if ( ( c->Chemical(2) / c->Area() ) < 0.1 )
    {
      c->EnlargeTargetArea( 0.02 * ( c->Area() ) );
    }
    else
    {
      c->EnlargeTargetArea( 0.2 * ( c->Area() ) );
    }
  }


  /// CAN BE REPLACED WITH ONE OF FOLLOWING MODULES... ///


  // LINEAR GROWTH WITH COUNTER: MODEL 1 (CF. TABLE S1 FOR PARAMETERS)

  if (mModel == 1 &&  c->CellType() != 0 )
  {
    time_now = m_mesh->getTime();

    if ( (((c->CellType()) >= 3) && ((c->CellType()) <= 7)) && (c->HasNeighborOfTypeZero() || ( ( c->NumberOfDivisions2() <= 2 ) && ( time_now - c->GetDivisionTime() ) >= 0)) )
    {
      c->EnlargeTargetArea( 2 );
    }
    else if ( (((c->CellType()) < 3) || ((c->CellType()) > 7)) && (c->HasNeighborOfTypeZero() || ( ( c->NumberOfDivisions2() <= 2 ) && ( time_now - c->GetDivisionTime() ) >= 0)) )
    {
      c->EnlargeTargetArea( 2 );
    }
    else if ( ( time_now - c->GetDivisionTime() ) <= 4320. && ( c->NumberOfDivisions2() > 2 ) )
    {
      c->EnlargeTargetArea( 20 );
    }
  }


  // EXPONENTIAL GROWTH WITH TIMER: MODEL 2,4 (CF. TABLE S1 FOR PARAMETERS)

  if ((mModel == 1 || mModel == 4) &&
      c->CellType() != 0 )
  {
    time_now = m_mesh->getTime();

    if ( (((c->CellType()) >= 3) && ((c->CellType()) <= 7)) && (c->HasNeighborOfTypeZero() || ( ( time_now - c->GetDivisionTime() ) <= 3240. && ( time_now - c->GetDivisionTime() ) >= 0)) )
    {
      c->EnlargeTargetArea( 0.02 * ( c->Area() ) );
    }
    else if ( (((c->CellType()) < 3) || ((c->CellType()) > 7)) && (c->HasNeighborOfTypeZero() || ( ( time_now - c->GetDivisionTime() ) <= 3240. && ( time_now - c->GetDivisionTime() ) >= 0)) )
    {
      c->EnlargeTargetArea( 0.02 * ( c->Area() ) );
    }
    else if ( ( time_now - c->GetDivisionTime() ) <= 4320. )
    {
      c->EnlargeTargetArea( 0.2 * ( c->Area() ) );
    }
  }


  // EXPONENTIAL GROWTH WITH COUNTER: MODEL 3,5,6,7 (CF. TABLE S1 FOR PARAMETERS)

  if ((mModel == 3 || mModel == 5 || mModel == 6 || mModel == 7) &&
      c->CellType() != 0 )
  {
    time_now = m_mesh->getTime();

    if ( (((c->CellType()) >= 3) && ((c->CellType()) <= 7)) && (c->HasNeighborOfTypeZero() || ( ( c->NumberOfDivisions2() <= 3 ) && ( time_now - c->GetDivisionTime() ) >= 0)) )
    {
      c->EnlargeTargetArea( 0.02 * ( c->Area() ) );
    }
    else if ( (((c->CellType()) < 3) || ((c->CellType()) > 7)) && (c->HasNeighborOfTypeZero() || ( ( c->NumberOfDivisions2() <= 3 ) && ( time_now - c->GetDivisionTime() ) >= 0)) )
    {
      c->EnlargeTargetArea( 0.02 * ( c->Area() ) );
    }
    else if ( ( time_now - c->GetDivisionTime() ) <= 4320. && ( c->NumberOfDivisions2() > 3 ) )
    {
      c->EnlargeTargetArea( 0.2 * ( c->Area() ) );
    }
  }


  // RULER + SIZER(CELL CYCLE): MODEL 8,9 (CF. TABLE S1 FOR PARAMETERS)

  if (mModel == 8 || mModel == 9)
  {
    time_now = m_mesh->getTime();
    double CCnoise2 = 1 /*+ (RANDOM()-0.5)/5*/;

    if ( ( tip_position - (c->Centroid().y) ) < /*180.*/240.)
    {
      if ( c->CellType() != 0 )
      {
        c->EnlargeTargetArea( 0.02 * ( c->Area() ) );
      }
    }
    else if ( ( tip_position - (c->Centroid().y) ) < 750. )
    {
      if ( c->CellType() != 0 )
      {
        c->EnlargeTargetArea( 0.2 * ( c->Area() ) );
      }
    }
  }

  // MODEL 10 (CF. TABLE S1 FOR PARAMETERS)
  if (mModel == 10)

  {
    time_now = m_mesh->getTime();

    if (time_now <= 4320)
    {
      if ( c->CellType() != 0 )
      {
        if ( (((c->CellType()) >= 3) && ((c->CellType()) <= 7)) && (c->HasNeighborOfTypeZero() || ( ( c->NumberOfDivisions2() <= 2 ) && ( time_now - c->GetDivisionTime() ) >= 0)) )
        {
          c->EnlargeTargetArea( 0.02 * ( c->Area() ) );
        }
        else if ( (((c->CellType()) < 3) || ((c->CellType()) > 7)) && (c->HasNeighborOfTypeZero() || ( ( c->NumberOfDivisions2() <= 2 ) && ( time_now - c->GetDivisionTime() ) >= 0)) )
        {
          c->EnlargeTargetArea( 0.02 * ( c->Area() ) );
        }
        else if ( ( time_now - c->GetDivisionTime() ) <= 44400. && ( c->NumberOfDivisions2() > 2 ) )
        {
          c->EnlargeTargetArea( 0.02 * ( c->Area() ) );
        }
      }
    }

    if ( time_now > 4320 )
    {
      if ( c->CellType() != 0 && ( c->NumberOfDivisions() > 2) )
      {
        if ( ( c->Chemical(0) / c->Area() ) > 13.5 )
        {
          c->EnlargeTargetArea( 0.02 * ( c->Area() ) );
        }
        else if ( ( c->Chemical(0) / c->Area() ) > 8.8 )
        {
          c->EnlargeTargetArea( 0.2 * ( c->Area() ) );
        }
      }
    }
  }

  // MODEL 11 (CF. TABLE S1 FOR PARAMETERS)
  if (mModel == 11)

  {
    time_now = m_mesh->getTime();

    if (time_now <= 4320)
    {
      if ( c->CellType() != 0 )
      {
        if ( (((c->CellType()) >= 3) && ((c->CellType()) <= 7)) && (c->HasNeighborOfTypeZero() || ( ( c->NumberOfDivisions2() <= 2 ) && ( time_now - c->GetDivisionTime() ) >= 0)) )
        {
          c->EnlargeTargetArea( 0.02 * ( c->Area() ) );
        }
        else if ( (((c->CellType()) < 3) || ((c->CellType()) > 7)) && (c->HasNeighborOfTypeZero() || ( ( c->NumberOfDivisions2() <= 2 ) && ( time_now - c->GetDivisionTime() ) >= 0)) )
        {
          c->EnlargeTargetArea( 0.02 * ( c->Area() ) );
        }
        else if ( ( time_now - c->GetDivisionTime() ) <= 44400. && ( c->NumberOfDivisions2() > 2 ) )
        {
          c->EnlargeTargetArea( 0.02 * ( c->Area() ) );
        }
      }
    }

    if ( time_now > 4320 )
    {
      if ( c->CellType() != 0 && (c->NumberOfDivisions() > 2) )
      {
        if ( c->CellType() == 3 || c->CellType() == 7 )
        {
          if ( ( c->Chemical(0) / c->Area() ) > 8. )
          {
            c->EnlargeTargetArea( 0.02 * ( c->Area() ) );
          }
          else if ( ( c->Chemical(0) / c->Area() ) > 6. )
          {
            c->EnlargeTargetArea( 0.2 * ( c->Area() ) );
          }
        }
        else
        {
          c->SetTargetArea(c->Area());//files passively following growth
        }
      }
    }
  }


  ///////////////////// **DIVISION RULES** /////////////////////////


  //  MODEL 12 (CF. TABLE S1 FOR PARAMETERS)
  if (mModel == 12)

  {
    double CCnoise2 = 1 + (RANDOM()-0.5)/5;

    if ( c->Index() > 60 )//avoids some artefacts from starting up
    {
      if ( (( c->Chemical(2) / c->Area() ) < 0.1) && (( c->Chemical(3) / c->Area() ) >= 0.7))
      {
        if ( ( c->CellType() > 2 && c->CellType() < 8 ) && ( c->Area() >= ( 200. * CCnoise2 ) ) && ( ( time_now - time_start ) >= 10. ) )
        {
          c->DivideOverAxis(Vector(1,0,0));
        }
        if ( c->CellType() != 0 && ( c->CellType() > 7 || c->CellType() < 3 ) && ( c->Area() >= ( 400.* CCnoise2 ) ) && ( ( time_now - time_start ) >= 10. ) )
        {
          c->DivideOverAxis(Vector(1,0,0));
        }
      }
    }
  }

  /// CAN BE REPLACED WITH ONE OF FOLLOWING MODULES... ///


  // COUNTER+TIMER(CELL CYCLE): MODELS 1, 3, 5 (CF. TABLE S1 FOR PARAMETERS)
  if (mModel == 1 || mModel == 3 || mModel == 5)

  {
    time_now = m_mesh->getTime();
    double CCnoise = 1 + (RANDOM()-0.5)/2;//to add noise at the start
    double CCnoise2 = 1 + (RANDOM()-0.5)/5;

    if ( (c->CellType() != 0) && (c->HasNeighborOfTypeZero()) && (time_now  == 0 || time_now == 30 ))
    {
      c->DivideOverAxis(Vector(1,0,0));
      c->SetDivCounter2( 0 );
    }

    if (( c->Index() >= 108 && c->Index() <= 119) && (c->NumberOfDivisions() == 2 ) && ( time_now - c->GetCellCycleTime() ) >= ( 1080. * CCnoise) )
    {
      c->SetDivCounter2( 0 );
      c->SetDivisionTime( m_mesh->getTime() );
      c->DivideOverAxis(Vector(1,0,0));
    }

    else
    {
      if ( ( c->CellType() > 3 && c->CellType() < 7 ) && ( ( time_now - c->GetCellCycleTime() ) >= ( 1080.* CCnoise2 ) ) )
      {
        if ( c->HasNeighborOfTypeZero() )
        {
          c->SetDivCounter2( 0 );
          c->SetDivisionTime( m_mesh->getTime() );
          c->DivideOverAxis(Vector(1,0,0));
        }
        else if ( ( time_now - c->GetDivisionTime() ) >= 0 && ( c->NumberOfDivisions2() <= 2 ) )
        {
          c->DivideOverAxis(Vector(1,0,0));
        }
      }
      if ( c->CellType() != 0 && ( c->CellType() > 6 || c->CellType() < 4 ) && ( ( time_now - c->GetCellCycleTime() ) >= ( 1080. * CCnoise2) ) )
      {
        if ( c->HasNeighborOfTypeZero() )
        {
          c->SetDivCounter2( 0 );
          c->SetDivisionTime( m_mesh->getTime() );
          c->DivideOverAxis(Vector(1,0,0));
        }
        else if ( time_now - c->GetDivisionTime() >= 0 && ( c->NumberOfDivisions2() <= 2) )
        {
          c->DivideOverAxis(Vector(1,0,0));
        }
      }
    }
  }

  // TIMER + TIMER (CELL CYCLE): MODELS 2, 4 (CF. TABLE S1 FOR PARAMETERS)
  if (mModel == 2 || mModel == 4)
  {
    time_now = m_mesh->getTime();
    double CCnoise = 1 + (RANDOM()-0.5)/2;
    double CCnoise2 = 1 + (RANDOM()-0.5)/5;//to add noise overall

    if ( (c->CellType() != 0) && (c->HasNeighborOfTypeZero()) && (time_now  == 0 || time_now == 30 ))
    {
      c->DivideOverAxis(Vector(1,0,0));
      c->SetDivCounter2( 0 );
    }

    if (( c->Index() >= 108 && c->Index() <= 119) && (c->NumberOfDivisions() == 2 ) && ( time_now - c->GetCellCycleTime() ) >= ( 1080. * CCnoise) )
    {
      c->SetDivCounter2( 0 );
      c->SetDivisionTime( m_mesh->getTime() );
      c->DivideOverAxis(Vector(1,0,0));
    }

    else
    {
      if ( ( c->CellType() > 3 && c->CellType() < 7 ) && ( time_now - c->GetCellCycleTime() >= ( 1080 * CCnoise2 ) ) )
      {
        if ( c->HasNeighborOfTypeZero() )
        {
          c->SetDivisionTime( m_mesh->getTime() );
          c->DivideOverAxis(Vector(1,0,0));
        }
        else if ( ( time_now - ( c->GetDivisionTime() ) ) >= 0 && ( ( time_now - c->GetDivisionTime() ) ) <= 3240 )
        {
          c->DivideOverAxis(Vector(1,0,0));
        }
      }
      if ( c->CellType() != 0 && ( c->CellType() > 6 || c->CellType() < 4 ) && ( time_now - c->GetCellCycleTime() >= ( 1080 * CCnoise2 ) ) )
      {
        if ( c->HasNeighborOfTypeZero() )
        {
          c->SetDivisionTime( m_mesh->getTime() );
          c->DivideOverAxis(Vector(1,0,0));
        }
        else if ( (time_now - c->GetDivisionTime() >= 0) && ( time_now - c->GetDivisionTime() ) <= 3240 )
        {
          c->DivideOverAxis(Vector(1,0,0));
        }
      }
    }
  }

  // COUNTER + SIZER(CELL CYCLE): MODELS 6,7 (CF. TABLE S1 FOR PARAMETERS)
  if (mModel == 6 || mModel == 7)
  {
    time_now = m_mesh->getTime();
    double CCnoise = 1 + (RANDOM()-0.5)/2;
    double CCnoise2 = 1 + (RANDOM()-0.5)/5;

    if ( (c->CellType() != 0) && (c->HasNeighborOfTypeZero()) && (time_now  == 0 || time_now == 30 ))
    {
      c->DivideOverAxis(Vector(1,0,0));
      c->SetDivCounter2( 0 );
    }

    if (( c->Index() >= 108 && c->Index() <= 119) && (c->NumberOfDivisions() == 2 ) )
    {
      if ( c->CellType() > 3 && c->CellType() < 7 && c->Area() >= 96. * CCnoise)
      {
        c->SetDivCounter2( 0 );
        c->SetDivisionTime( m_mesh->getTime() );
        c->DivideOverAxis(Vector(1,0,0));
      }
      else if ( ( c->CellType() < 3 || c->CellType() > 7 ) && c->Area() >= 160. * CCnoise)
      {
        c->SetDivCounter2( 0 );
        c->SetDivisionTime( m_mesh->getTime() );
        c->DivideOverAxis(Vector(1,0,0));
      }
      else if ( ( c->CellType() == 3 || c->CellType() == 7 ) && c->Area() >= 80. * CCnoise)
      {
        c->SetDivCounter2( 0 );
        c->SetDivisionTime( m_mesh->getTime() );
        c->DivideOverAxis(Vector(1,0,0));
      }
    }

    else
    {
      if ( ( c->CellType() > 3 && c->CellType() < 7 ) && ( c->Area() >= ( 96. * CCnoise2) ) )
      {
        if ( c->HasNeighborOfTypeZero() )
        {
          c->SetDivCounter2( 0 );
          c->SetDivisionTime( m_mesh->getTime() );
          c->DivideOverAxis(Vector(1,0,0));
        }
        else if ( ( time_now - c->GetDivisionTime() ) >= 0 && ( c->NumberOfDivisions2() <= 3 ) )
        {
          c->DivideOverAxis(Vector(1,0,0));
        }
      }
      if ( c->CellType() != 0 && ( c->CellType() > 7 || c->CellType() < 3 ) && ( c->Area() >= ( 160.* CCnoise2) ) )
      {
        if ( c->HasNeighborOfTypeZero() )
        {
          c->SetDivCounter2( 0 );
          c->SetDivisionTime( m_mesh->getTime() );
          c->DivideOverAxis(Vector(1,0,0));
        }
        else if ( time_now - c->GetDivisionTime() >= 0 && ( c->NumberOfDivisions2() <= 3) )
        {
          c->DivideOverAxis(Vector(1,0,0));
        }
      }
      if ( ( c->CellType() == 7 || c->CellType() == 3 ) && ( c->Area() >= ( 80.* CCnoise2 ) ) )
      {
        if ( c->HasNeighborOfTypeZero() )
        {
          c->SetDivCounter2( 0 );
          c->SetDivisionTime( m_mesh->getTime() );
          c->DivideOverAxis(Vector(1,0,0));
        }
        else if ( time_now - c->GetDivisionTime() >= 0 && ( c->NumberOfDivisions2() <= 3) )
        {
          c->DivideOverAxis(Vector(1,0,0));
        }
      }
    }
  }

  // RULER + SIZER(CELL CYCLE): MODEL 8 (CF. TABLE S1 FOR PARAMETERS)
  if (mModel == 8)
  {
    time_now = m_mesh->getTime();
    double CCnoise2 = 1 + (RANDOM()-0.5)/5;

    if ( ( tip_position - (c->Centroid().y) ) < 180.)
    {
      if ( ( c->CellType() > 3 && c->CellType() < 7 ) && ( c->Area() >= ( 120. * CCnoise2) ) )
      {
        c->DivideOverAxis(Vector(1,0,0));
      }
      if ( c->CellType() != 0 && ( c->CellType() > 7 || c->CellType() < 3 ) && ( c->Area() >= ( 240.* CCnoise2) ) )
      {
        c->DivideOverAxis(Vector(1,0,0));
      }
      if ( ( c->CellType() == 7 || c->CellType() == 3 ) && ( c->Area() >= ( 144. * CCnoise2 ) ) )
      {
        c->DivideOverAxis(Vector(1,0,0));
      }
    }
  }

  // RULER + SIZER(CELL CYCLE): MODEL 9 (CF. TABLE S1 FOR PARAMETERS)
  if (mModel == 9)

  {
    time_now = m_mesh->getTime();
    double CCnoise2 = 1 + (RANDOM()-0.5)/5;

    if ( ( tip_position - (c->Centroid().y) ) < 240.)
    {
      if ( ( c->CellType() != 0 ) && ( c->Area() >= ( 400. * CCnoise2) ) )
      {
        c->DivideOverAxis(Vector(1,0,0));
      }
    }

  }

  // MODEL 10, 11 (CF. TABLE S1 FOR PARAMETERS)
  if (mModel == 10 || mModel == 11)
  {
    time_now = m_mesh->getTime();

    if (time_now <= 4320)
    {
      double CCnoise = 1 + (RANDOM()-0.5)/2;
      double CCnoise2 = 1 /*+ (RANDOM()-0.5)/5*/;

      if ( (c->CellType() != 0) && (c->HasNeighborOfTypeZero()) && (time_now  == 0 || time_now == 30 ))
      {
        c->DivideOverAxis(Vector(1,0,0));
        c->SetDivCounter2( 0 );
      }

      if (( c->Index() >= 108 && c->Index() <= 119) && (c->NumberOfDivisions() == 2 ) )
      {
        if ( c->CellType() > 3 && c->CellType() < 7 && c->Area() >= 96. * CCnoise)
        {
          c->SetDivCounter2( 0 );
          c->SetDivisionTime( m_mesh->getTime() );
          c->DivideOverAxis(Vector(1,0,0));
        }
        else if ( ( c->CellType() < 3 || c->CellType() > 7 ) && c->Area() >= 160. * CCnoise)
        {
          c->SetDivCounter2( 0 );
          c->SetDivisionTime( m_mesh->getTime() );
          c->DivideOverAxis(Vector(1,0,0));
        }
        else if ( ( c->CellType() == 3 || c->CellType() == 7 ) && c->Area() >= 80. * CCnoise)
        {
          c->SetDivCounter2( 0 );
          c->SetDivisionTime( m_mesh->getTime() );
          c->DivideOverAxis(Vector(1,0,0));
        }
      }

      else
      {
        if ( ( c->CellType() > 3 && c->CellType() < 7 ) && ( c->Area() >= ( 96. * CCnoise2) ) )
        {
          if ( c->HasNeighborOfTypeZero() )
          {
            c->SetDivCounter2( 0 );
            c->SetDivisionTime( m_mesh->getTime() );
            c->DivideOverAxis(Vector(1,0,0));
          }
          else if ( ( time_now - c->GetDivisionTime() ) >= 0 && ( c->NumberOfDivisions2() <= 2) )
          {
            c->DivideOverAxis(Vector(1,0,0));
          }
        }
        if ( c->CellType() != 0 && ( c->CellType() > 7 || c->CellType() < 3 ) && ( c->Area() >= ( 160.* CCnoise2) ) )
        {
          if ( c->HasNeighborOfTypeZero() )
          {
            c->SetDivCounter2( 0 );
            c->SetDivisionTime( m_mesh->getTime() );
            c->DivideOverAxis(Vector(1,0,0));
          }
          else if ( time_now - c->GetDivisionTime() >= 0 && ( c->NumberOfDivisions2() <= 2) )
          {
            c->DivideOverAxis(Vector(1,0,0));
          }
        }
        if ( ( c->CellType() == 7 || c->CellType() == 3 ) && ( c->Area() >= ( 80.* CCnoise2) ) )
        {
          if ( c->HasNeighborOfTypeZero() )
          {
            c->SetDivCounter2( 0 );
            c->SetDivisionTime( m_mesh->getTime() );
            c->DivideOverAxis(Vector(1,0,0));
          }
          else if ( time_now - c->GetDivisionTime() >= 0 && ( c->NumberOfDivisions2() <= 2) )
          {
            c->DivideOverAxis(Vector(1,0,0));
          }
        }
      }
    }

    if ( time_now > 4320 )
    {
      double CCnoise2 = 1 /*+ (RANDOM()-0.5)/5*/;

      //if ( ( c->Chemical(0) / c->Area() ) > 13.5 )// MODEL 10
      if ( ( tip_position - (c->Centroid().y) ) < 180.)// MODEL 11
      {
        if ( ( c->CellType() > 3 && c->CellType() < 7 ) && ( c->Area() >= ( 96. * CCnoise2 ) ) )
        {
          c->DivideOverAxis(Vector(1,0,0));
        }
        if ( c->CellType() != 0 && ( c->CellType() > 7 || c->CellType() < 3 ) && ( c->Area() >= ( 160.* CCnoise2 ) ) )
        {
          c->DivideOverAxis(Vector(1,0,0));
        }
        if ( ( c->CellType() == 7 || c->CellType() == 3 ) && ( c->Area() >= ( 80.* CCnoise2 ) ) )
        {
          c->DivideOverAxis(Vector(1,0,0));
        }
      }
    }

  }
}

void Wortel::CelltoCellTransport(Wall* w, double* dchem_c1, double* dchem_c2)
{

  //Parameter&      par = Parameter::getInstance();

  if ( (w->C1()->BoundaryPolP()) || (w->C2()->BoundaryPolP()) )
  {
    return;
  }

  double const apoplast_thickness = par->apoplast_thickness;
  double const phi = (w->Length() / apoplast_thickness) * par->D[0] * ( ( w->C2()->Chemical(0) / (w->C2()->Area()) )
                                                                        - ( w->C1()->Chemical(0) / (w->C1()->Area()) ) );

  dchem_c1[0] += phi;
  dchem_c2[0] -= phi;

  // Active fluxes (PIN1 mediated transport)
  double const k_import = par->k_import;
  double const k_export = par->transport;
  double const kmshy = 0.1;

  // transport from cell 1 to cell 2
  double const shy12 = kmshy  / ( kmshy + ( w->C1()->Chemical(2) ) / (w->C1()->Area()) );
  double const trans12 = w->Length() * ( w->C1()->Chemical(0) / (w->C1()->Area()) ) * (k_export * w->Transporters1(1) * shy12 + k_import);

  // transport from cell 2 to cell 1
  double const shy21 = kmshy  / ( kmshy + ( w->C2()->Chemical(2) ) / (w->C2()->Area()) );
  double const trans21 = w->Length() * ( w->C2()->Chemical(0) / (w->C2()->Area()) ) * (k_export * w->Transporters2(1) * shy21 + k_import);

  dchem_c1[0] += (trans21 - trans12);
  dchem_c2[0] += (trans12 - trans21);


  //Second diffusive hormone - chemical '1'
  double phi2 = (w->Length() / apoplast_thickness) * 12 * ( ( w->C2()->Chemical(1) / (w->C2()->Area()) )//D[1]=12
                                                            - ( w->C1()->Chemical(1) / (w->C1()->Area()) ) );

  dchem_c1[1] += phi2 ;
  dchem_c2[1] -= phi2 ;

}


double Wortel::PINflux(CellBase* this_cell, CellBase* adjacent_cell, Wall* w)
{
  // PIN1 localization at wall
  // Note: chemical 0 is Auxin (intracellular storage only)
  //  PIN1 is Chemical 1 (both in walls and intracellular storage)
  //! \f$ \frac{d Pij/dt}{dt} = k_1 A_j \frac{P_i}{L_ij} - k_2 P_{ij} \f$
  // Note that Pij is measured in term of concentration (mol/L)
  // Pi in terms of quantity (mol)

  // Equations as in Merks et al., Trends in Plant Science 2007

  // calculate PIN translocation rate from cell to membrane

  double const adj_auxin = adjacent_cell->Chemical(0);
  double const receptor_level = adj_auxin * par->r / (par->kr + adj_auxin);
  double pin_atwall; // pick the correct side of the Wall

  // note: pin_flux is net flux from endosome to wall
  double const pin_flux = 0;/*par->k1 * this_cell->Chemical(1) * receptor_level / (par->km + this_cell->Chemical(1))
                      - par->k2 * pin_atwall;*/

  return pin_flux;
}


void Wortel::WallDynamics(Wall* w, double* dw1, double* dw2)
{
  // add biochemical networks for reactions occurring at walls here

}

double Wortel::Michaelis(double M1, double J1, double K1, double S1)
{
  double rate = 0;
  return rate = M1 * S1 / ( J1 + K1 );

}

double Wortel::Goldbeter(double A1, double A2, double A3, double A4)
{
  double rate = 0;
  double BB = A2 - A1 + A3 * A2 + A4 * A1 ;
  return rate = 2.0 * A4 * A1 / ( BB + pow ( ( pow ( BB, 2.0 ) - 4.0 * ( A2 - A1 ) * A4 * A1 ), .5 ) ) ;
}

double Wortel::Hill(double Vm, double Km, double h, double S)
{
  double rate = 0;
  return rate = Vm * pow ( S, h ) / ( pow ( S, h ) + pow ( Km, h ) ) ;
}

void Wortel::CellDynamics(CellBase* c, double* dchem)
{
  // add biochemical networks for intracellular reactions here



  ////auxin and cytokinin dynamics (MODEL 12)
  if (mModel != 12) return;

  double kmauxCK = 100.;
  if (c->Index() >= 2 && c->Index() <= 9 )
  {
    dchem[0] = 100000 + (c->Area()) * par->aux1prod - par->aux_breakdown * (c->Chemical(0)) ;
    dchem[1] = 25 + ( 0.01 * (c->Area()) * kmauxCK / ( kmauxCK + (c->Chemical(0) / (c->Area())) ) ) - 0.0005 * (c->Chemical(1)) ;
  }
  else if ( c->Index() == 0 || c->Index() <= 1 || c->Index() >= 10 || c->Index() <= 11 )
  {
    dchem[0] = -100000 + (c->Area()) * par->aux1prod - par->aux_breakdown * (c->Chemical(0)) ;
    dchem[1] = -25 + ( 0.01 * (c->Area()) * kmauxCK / ( kmauxCK + (c->Chemical(0) / (c->Area())) ) ) - 0.0005 * (c->Chemical(1)) ;
  }
  else
  {
    dchem[0] = (c->Area()) * par->aux1prod - par->aux_breakdown * (c->Chemical(0)) ;
    dchem[1] = ( 0.01 * (c->Area()) * kmauxCK / ( kmauxCK + (c->Chemical(0) / (c->Area())) ) ) - 0.0005 * (c->Chemical(1)) ;
  }

  ////SHY2 dynamics
  double kmaux = 100;
  if( ( c->Chemical(1) / (c->Area()) )  > 2.5 )
  {
    dchem[2] = 1 * (c->Area()) - (c->Chemical(2)) * ( 0.001 + 0.1 * ( (c->Chemical(0) / c->Area()) / (kmaux + c->Chemical(0) / c->Area()) ) );
  }
  else
  {
    dchem[2] = - 0.001 * (c->Chemical(2));
  }

  ////GA dynamics

  if ( c->CellType() > 2 && c->CellType() < 8 )
  {
    dchem[3] = 100 - 0.1 * c->Chemical(3);
  }
  else
  {
    dchem[3] = 200 - 0.1 * c->Chemical(3);//wider cells have proportionally more GA production
  }

}

Q_EXPORT_PLUGIN2( wortel, Wortel )
