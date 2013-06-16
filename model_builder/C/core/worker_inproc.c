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

void *worker_routine_smc_inproc(void *params) 
{
    int j, n, np1, t0, t1;
    int id;
    char str[STR_BUFFSIZE];

    struct s_thread_smc *p = (struct s_thread_smc *) params;

    struct s_data *p_data = p->p_data;
    struct s_par *p_par = p->p_par;
    struct s_calc *p_calc = p->p_calc;
    struct s_X ***D_J_p_X = p->D_J_p_X;
    struct s_likelihood *p_like = p->p_like;

    plom_f_pred_t f_pred = get_f_pred(p_data->implementation, p_data->noises_off);

    // Socket to server controller
    void *controller = zmq_socket (p->context, ZMQ_SUB);
    zmq_connect (controller, "inproc://server_controller");
    zmq_setsockopt (controller, ZMQ_SUBSCRIBE, "", 0);

    //  Socket to receive messages (particles) from the server
    void *receiver = zmq_socket (p->context, ZMQ_PULL);   
    zmq_connect (receiver, "inproc://server_sender");

    //  Socket to send messages (results) to the server
    void *sender = zmq_socket (p->context, ZMQ_PUSH);
    zmq_connect (sender, "inproc://server_receiver");

    // ready !    
    zmq_send(sender, &p->thread_id, sizeof (int), 0);

    zmq_pollitem_t items [] = {
        { receiver, 0, ZMQ_POLLIN, 0 },
        { controller, 0, ZMQ_POLLIN, 0 }
    };

    while (1) {
        zmq_poll (items, 2, -1);
        if (items [0].revents & ZMQ_POLLIN) {

	    zmq_recv(receiver, &id, sizeof (int), 0);

	    n = p_calc->current_n;
	    np1 = n+1;
	    t0 = p_data->times[n];
	    t1 = p_data->times[np1];

	    int J_start = id * p->J_chunk;
	    int J_end = (id+1 == p_calc->n_threads) ? p->J : (id+1)*p->J_chunk;	  
	   
	    for(j=J_start; j<J_end; j++ ){
		reset_inc(D_J_p_X[np1][j], p_data);
		(*f_pred)(D_J_p_X[np1][j], t0, t1, p_par, p_data, p_calc);

		proj2obs(D_J_p_X[np1][j], p_data);

		if(p_data->data_ind[n]->n_nonan) {
		    p_like->weights[j] = exp(get_log_likelihood(D_J_p_X[np1][j], p_par, p_data, p_calc));
		}
	    }

	    //send back id of the batch of particles now integrated
	    zmq_send(sender, &p->thread_id, sizeof (int), 0);
	}

        //controller commands:
        if (items [1].revents & ZMQ_POLLIN) {
	    char buf [256];
	    zmq_recv(controller, buf, 256, 0);

	    snprintf(str, STR_BUFFSIZE, "worker %d: controller sent: %s", p->thread_id, buf);
	    print_log(str);

            if(strcmp(buf, "KILL") == 0) {
                break;  //  Exit loop
            }
        }
    }

    zmq_close (receiver);
    zmq_close (sender);
    zmq_close (controller);

    snprintf(str, STR_BUFFSIZE, "thread %d done", p->thread_id);
    print_log(str);

    return NULL;
}


void *worker_routine_mif_inproc(void *params) 
{
    int j, n, np1, t0, t1;
    int id;
    char str[STR_BUFFSIZE];

    struct s_thread_mif *p = (struct s_thread_mif *) params;

    struct s_data *p_data = p->p_data;
    struct s_par **J_p_par = p->J_p_par;
    struct s_calc *p_calc = p->p_calc;
    struct s_X ***J_p_X = p->J_p_X;
    struct s_likelihood *p_like = p->p_like;

    plom_f_pred_t f_pred = get_f_pred(p_data->implementation, p_data->noises_off);


    void *controller = zmq_socket (p->context, ZMQ_SUB);
    zmq_connect (controller, "inproc://server_controller");
    zmq_setsockopt (controller, ZMQ_SUBSCRIBE, "", 0);

    void *receiver = zmq_socket (p->context, ZMQ_PULL);   
    zmq_connect (receiver, "inproc://server_sender");

    void *sender = zmq_socket (p->context, ZMQ_PUSH);
    zmq_connect (sender, "inproc://server_receiver");

    // ready !    
    zmq_send(sender, &p->thread_id, sizeof (int), 0);


    zmq_pollitem_t items [] = {
        { receiver, 0, ZMQ_POLLIN, 0 },
        { controller, 0, ZMQ_POLLIN, 0 }
    };

    while (1) {
        zmq_poll (items, 2, -1);
        if (items [0].revents & ZMQ_POLLIN) {
	   
	    zmq_recv(receiver, &id, sizeof (int), 0);

	    n = p_calc->current_n;
	    np1 = n+1;
	    t0 = p_data->times[n];
	    t1 = p_data->times[np1];

	    int J_start = id * p->J_chunk;
	    int J_end = (id+1 == p_calc->n_threads) ? p->J : (id+1)*p->J_chunk;	  

	    for(j=J_start; j<J_end; j++ ){
		reset_inc((*J_p_X)[j], p_data);

		f_pred((*J_p_X)[j], t0, t1, J_p_par[j], p_data, p_calc);

		if(p_data->data_ind[n]->n_nonan) {
		    proj2obs((*J_p_X)[j], p_data);
		    p_like->weights[j] = exp(get_log_likelihood((*J_p_X)[j], J_p_par[j], p_data, p_calc));
		}
	    }

	    //send back id of the batch of particles now integrated
	    zmq_send(sender, &p->thread_id, sizeof (int), 0);
	}

        //controller commands:
        if (items [1].revents & ZMQ_POLLIN) {
	    char buf [256];
	    zmq_recv(controller, buf, 256, 0);           

	    snprintf(str, STR_BUFFSIZE, "worker %d: controller sent: %s", p->thread_id, buf);
	    print_log(str);

            if(strcmp(buf, "KILL") == 0) {
                break;  //  Exit loop
            }
        }
    }

    zmq_close (receiver);
    zmq_close (sender);
    zmq_close (controller);

    snprintf(str, STR_BUFFSIZE, "thread %d done", p->thread_id);
    print_log(str);

    return NULL;
}


void *worker_routine_predict_inproc(void *params) 
{
    int j;
    int k, kp1, id;
    char str[STR_BUFFSIZE];

    struct s_thread_predict *p = (struct s_thread_predict *) params;

    struct s_data *p_data = p->p_data;
    struct s_par **J_p_par = p->J_p_par;
    struct s_calc *p_calc = p->p_calc;
    struct s_X **J_p_X = p->J_p_X;   

    plom_f_pred_t f_pred = get_f_pred(p_data->implementation, p_data->noises_off);

    void *controller = zmq_socket (p->context, ZMQ_SUB);
    zmq_connect (controller, "inproc://server_controller");
    zmq_setsockopt (controller, ZMQ_SUBSCRIBE, "", 0);

    void *receiver = zmq_socket (p->context, ZMQ_PULL);   
    zmq_connect (receiver, "inproc://server_sender");

    void *sender = zmq_socket (p->context, ZMQ_PUSH);
    zmq_connect (sender, "inproc://server_receiver");

    // ready !    
    zmq_send(sender, &p->thread_id, sizeof (int), 0);

    zmq_pollitem_t items [] = {
        { receiver, 0, ZMQ_POLLIN, 0 },
        { controller, 0, ZMQ_POLLIN, 0 }
    };

    while (1) {
        zmq_poll (items, 2, -1);
        if (items [0].revents & ZMQ_POLLIN) {

	    zmq_recv(receiver, &id, sizeof (int), 0);
	    zmq_recv(receiver, &k, sizeof (int), 0);
	    kp1 = k+1;

	    int J_start = id * p->J_chunk;
	    int J_end = (id+1 == p_calc->n_threads) ? p->J : (id+1)*p->J_chunk;	  

	    for(j=J_start; j<J_end; j++ ){
		reset_inc(J_p_X[j], p_data);
		f_pred(J_p_X[j], k, kp1, J_p_par[j], p_data, p_calc);
		proj2obs(J_p_X[j], p_data);
	    }

	    //send back id of the batch of particles now integrated
	    zmq_send(sender, &p->thread_id, sizeof (int), 0);
	}

        //controller commands:
        if (items [1].revents & ZMQ_POLLIN) {
	    char buf [256];
	    zmq_recv(controller, buf, 256, 0);           	    

	    snprintf(str, STR_BUFFSIZE, "worker %d: controller sent: %s", p->thread_id, buf);
	    print_log(str);

            if(strcmp(buf, "KILL") == 0) {
                break;  //  Exit loop
            }
        }
    }

    zmq_close (receiver);
    zmq_close (sender);
    zmq_close (controller);

    snprintf(str, STR_BUFFSIZE, "thread %d done", p->thread_id);
    print_log(str);

    return NULL;
}
