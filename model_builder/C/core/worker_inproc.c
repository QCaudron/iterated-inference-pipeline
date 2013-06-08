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

void *worker_routine_smc_inproc (void *params) 
{
    int j, nn, nnp1, t1;
    int id;

    struct s_thread_smc *p = (struct s_thread_smc *) params;

    // Socket to server controller
    void *server_controller = zmq_socket (p->context, ZMQ_SUB);
    zmq_connect (server_controller, "inproc://plom_controller");
    zmq_setsockopt (server_controller, ZMQ_SUBSCRIBE, "", 0);

    //  Socket to receive messages (particles) from the server
    void *server_receiver = zmq_socket (p->context, ZMQ_PULL);   
    zmq_connect (server_receiver, "inproc://server_sender");

    //  Socket to send messages (results) to the server
    void *server_sender = zmq_socket (p->context, ZMQ_PUSH);
    zmq_connect (server_sender, "inproc://server_receiver");

    struct s_data *p_data = p->p_data;
    struct s_par *p_par = p->p_par;
    struct s_calc **calc = p->calc;
    struct s_X ***D_J_p_X = p->D_J_p_X;
    struct s_likelihood *p_like = p->p_like;
    int size_J = p->size_J;

    plom_f_pred_t f_pred = get_f_pred(p_data->implementation, p_data->noises_off);

    zmq_pollitem_t items [] = {
        { server_receiver, 0, ZMQ_POLLIN, 0 },
        { server_controller, 0, ZMQ_POLLIN, 0 }
    };

    while (1) {
        zmq_poll (items, 2, -1);
        if (items [0].revents & ZMQ_POLLIN) {

	    //particule to integrate
	    zmq_recv(server_receiver, &id, sizeof (int), 0);
	    //printf("rcv %d from thread %d\n", j, id);

	    for(j=id*size_J; j<(id+1)*size_J; j++ ){
		nn = calc[id]->current_nn;
		nnp1 = nn+1;
		t1 = p_data->times[calc[id]->current_n];

		reset_inc(D_J_p_X[nnp1][j], p_data);
		(*f_pred)(D_J_p_X[nnp1][j], nn, nnp1, p_par, p_data, calc[id]);

		proj2obs(D_J_p_X[nnp1][j], p_data);

		if(nnp1 == t1) {
		    p_like->weights[j] = exp(get_log_likelihood(D_J_p_X[nnp1][j], p_par, p_data, calc[id]));
		}
	    }

	    //send back paricle index (now integrated)
	    zmq_send(server_sender, &j, sizeof (int), 0);
	}

        //controller commands:
        if (items [1].revents & ZMQ_POLLIN) {
	    char buf [256];
	    zmq_recv(server_controller, buf, 256, 0);           
            if(strcmp(buf, "KILL") == 0) {
                printf("worker %d: controller sent: %s\n", p->thread_id, buf);
                break;  //  Exit loop
            }
        }
    }

    zmq_close (server_receiver);
    zmq_close (server_sender);
    zmq_close (server_controller);

    printf("thread %d done\n", p->thread_id);

    return NULL;
}
