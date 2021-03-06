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

void ask_update()
{
    json_t *root;
    root = json_pack("{s,s}", "flag", "upd");
    json_dumpf(root, stderr, 0); fprintf(stderr, "\n");
    fflush(stderr);
    json_decref(root);
}


/**
 * webApp only.
 * Block, waiting for an empty msg from stdin. This is used to prevent
 * the webserver and the client to be saturated by output comming from
 * the C code.
 */
void block()
{
    //TO DO: add a new flag in ask_update != upd for simple blocking.
    ask_update();
    char tmp[4];
    fgets(tmp, 4, stdin); //we will get "{}\n" when the websocket server is ready...
}


/**
 *   live update of the walk rates integrating data from the webApp.
 */
void update_walk_rates(struct s_best *p_best, struct s_data *p_data)
{
    ask_update();

    fflush(stderr); //just to be sure that ask_update has delivered its message

    json_t *theta = load_json();
//	print_json_on_stdout(root);

    if(json_object_size(theta)) { 
        load_best(p_best, p_data, theta, 0);
        print_err("data integrated");
    }

    json_decref(theta);
}
