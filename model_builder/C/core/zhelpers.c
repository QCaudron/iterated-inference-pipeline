/**************************************************************************
 *    This file is part of plom.
 *
 *    plom is free software: you can redistribute it and/or modify it
 *    under the terms of the GNU General Public License as published
 *    by the Free Software Foundation, either version 3 of the
 *    License, or (at your option) any later version.
 *
 *    plom is distributed in the hope that it will be useful, but
 *    WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public
 *    License along with plom.  If not, see
 *    <http://www.gnu.org/licenses/>.
 *************************************************************************/

#include "plom.h"


void send_par(void *socket, const struct s_par *p_par, struct s_data *p_data, int zmq_options)
{   
    int i;

    //send "natural" component of p_par
    for(i=0; i< p_par->size_natural; i++) {
	zmq_send(socket, p_par->natural[i], p_data->routers[i]->n_gp *sizeof (double), zmq_options);
    }

}


void recv_par(struct s_par *p_par, struct s_data *p_data, void *socket)
{
    int i;

    //natural
    for(i=0; i< p_par->size_natural; i++) {
	zmq_recv(socket, p_par->natural[i], p_data->routers[i]->n_gp * sizeof (double), 0);               
    }
}


void send_X(void *socket, const struct s_X *p_X, struct s_data *p_data, int zmq_options)
{
    int size_proj = N_PAR_SV*N_CAC + p_data->p_it_only_drift->nbtot + N_TS_INC_UNIQUE;

    //dt
    zmq_send(socket, &(p_X->dt), sizeof (double), ZMQ_SNDMORE);    

    //send obs
    zmq_send(socket, p_X->obs, N_TS * sizeof (double), ZMQ_SNDMORE);        
   
    //send proj
    zmq_send(socket, p_X->proj, size_proj * sizeof (double), zmq_options);            

}


void recv_X(struct s_X *p_X, struct s_data *p_data, void *socket)
{

    int size_proj = N_PAR_SV*N_CAC + p_data->p_it_only_drift->nbtot + N_TS_INC_UNIQUE;

    //dt
    zmq_recv(socket, &(p_X->dt), sizeof (double), 0);

    //obs
    zmq_recv(socket, p_X->obs, N_TS * sizeof (double), 0);        

    //proj
    zmq_recv(socket, p_X->proj, size_proj * sizeof (double), 0);        

}
