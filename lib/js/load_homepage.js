function load_viz_data(inst_name){

  // var inst_name = 'mult_view';
  var tmp_num;
  var cat_colors;
  cgm = {};

  resize_container();

  function make_clust(){


    var clust_name = 'homepage_phos.json'

    d3.json('json/'+clust_name, function(network_data){
      var args = $.extend(true, {}, default_args);
      args.root = '#container-id-1';
      args.network_data = network_data;
      args.row_tip_callback = gene_info

      cgm['phos'] = Clustergrammer(args);
      d3.select(cgm['phos'].params.root+' .wait_message').remove();
      cat_colors = cgm['phos'].params.viz.cat_colors;

      // enr_obj = Enrichr_request(cgm.clust);
      // enr_obj.enrichr_icon();

      $.unblockUI();

    });

    var clust_name = 'homepage_exp.json'

    d3.json('json/'+clust_name, function(network_data){
      var args = $.extend(true, {}, default_args);
      args.root = '#container-id-2';
      args.network_data = network_data;
      args.row_tip_callback = gene_info

      cgm['exp'] = Clustergrammer(args);
      d3.select(cgm['exp'].params.root+' .wait_message').remove();
      cat_colors = cgm['exp'].params.viz.cat_colors;

      enr_obj = Enrichr_request(cgm.exp);
      enr_obj.enrichr_icon();

      $.unblockUI();

    });

  }

  // make wait sign
  $.blockUI({ css: {
      border: 'none',
      padding: '15px',
      backgroundColor: '#000',
      '-webkit-border-radius': '10px',
      '-moz-border-radius': '10px',
      opacity: .8,
      color: '#fff'
  } });

  d3.select('.blockMsg').select('h1').text('Please wait...');

  var viz_size = {'width':1140, 'height':750};

  // define arguments object
  var default_args = {
    'show_tile_tooltips':true,
    'row_search_placeholder':'Gene',
    'sidebar_icons': false,
    'col_label_scale':1.5
  };

  $(document).ready(function(){
      $(this).scrollTop(0);
  });

  make_clust()

  d3.select(window).on('resize',function(){
    resize_container();

    _.each(cgm, function(inst_cgm){
      inst_cgm.resize_viz();
    })

  });

  window.onscroll = function() {

    var show_col_sim = 200;
    var show_row_sim = 1200;
    var hide_clust = 900;
    var hide_col_sim = 1800;
    var inst_scroll = $(document).scrollTop();

    // hide clust
    if (inst_scroll > hide_clust){
      d3.select('#container-id-1 .viz_svg')
        .style('display', 'none');
    } else {
      d3.select('#container-id-1 .viz_svg')
        .style('display', 'block');
    }

    // hide col sim mat
    if (inst_scroll > hide_col_sim || inst_scroll < show_col_sim){
      d3.select('#container-id-2 .viz_svg')
        .style('display', 'none');
    } else {
      d3.select('#container-id-2 .viz_svg')
        .style('display', 'block');
    }

  }



  function unblock(){
    $.unblockUI();
    // d3.selectAll('.row_label_group text').style('font-family','"Courier new"')
    // d3.selectAll('.col_label_text').style('font-family','"Courier new"')
  }

  function resize_container(){

    var container_width = d3.select('#wrap').style('width').replace('px','');

    var container_width = Number(container_width) - 30;

    d3.selectAll('.clustergrammer_container')
      .style('width', container_width+'px');

    // var container_height = 700;
    // d3.selectAll('.sim_mat_container')
    //   .style('height', container_height+'px');

  }

}


// save gene data
gene_data = {};

function gene_info(gene_symbol){

  if (_.has(gene_data, gene_symbol)){
    var inst_data = gene_data[gene_symbol];
    set_tooltip(inst_data)
  } else{
    setTimeout(hzome_get_request, 500, gene_symbol);
  }

  function hzome_get_request(gene_symbol){

    if ( d3.select(' .row_tip').classed(gene_symbol) ){

      gene_symbol = gene_symbol.split(' ')[0];
      gene_symbol = gene_symbol.split('-')[0];

      // var base_url = 'http://localhost:9000/clustergrammer/gene_info/';
      var base_url = 'http://amp.pharm.mssm.edu/Harmonizome/api/1.0/gene/'
      // var base_url = 'http://yan.1425mad.mssm.edu:31331/Harmonizome/api/1.0/gene/'

      console.log('hzome')
      var url = base_url + gene_symbol;

      $.get(url, function(data) {
        data = JSON.parse(data);

        // save data for repeated use
        gene_data[gene_symbol] = {}
        gene_data[gene_symbol].name = data.name;
        gene_data[gene_symbol].description = data.description;

        set_tooltip(data);

      });

    }

  }

  function set_tooltip(data){

    if (data.name != undefined){
      d3.select('.row_tip')
        .html(function(){
            var sym_name = gene_symbol + ': ' + data.name;
            var full_html = '<p>' + sym_name + '</p>' +  '<p>' +
              data.description + '</p>';
            return full_html;
        });
    }
  }

}





function update_viz_callback(enr_obj){

  cgm.clust.update_cats(enr_obj.cat_data);

  d3.select(cgm.clust.params.root+' .enr_title').remove();

  var enr_title = d3.select(cgm.clust.params.root+' .viz_svg')
    .append('g')
    .classed('enr_title', true)
    .attr('transform', function(){

      var trans = d3.select('.row_cat_label_container')
                    .attr('transform').split('(')[1].split(')')[0];
      x_offset = Number(trans.split(',')[0]) - 10;

      return 'translate('+ String(x_offset)+', 0)';

    });

  enr_title
    .append('rect')
    .attr('width', cgm.clust.params.viz.cat_room.row)
    .attr('height', 25)
    .attr('fill', 'white');

  var library_string = enr_obj.library.substring(0,40);
  enr_title
    .append('text')
    .attr('transform', 'translate(0, 17)')
    .text(library_string.replace(/_/g, ' '))
    .style('font-size', '15px')
    .attr('font-family', '"Helvetica Neue", Helvetica, Arial, sans-serif');

}

function make_enr_wait_circle(){
  var pos_x = 71;
  var pos_y = 25;

   var click_circle = d3.select(cgm.clust.params.root+' .viz_svg')
      .append('circle')
      .classed('enr_wait_circle', true)
      .attr('cx',pos_x)
      .attr('cy',pos_y)
      .attr('r',22)
      .style('stroke','#666666')
      .style('stroke-width','3px')
      .style('fill','white')
      .style('fill-opacity',0)
      .style('opacity', 0);
}


function animate_wait() {
  var repeat_time = 700;
  var max_opacity = 0.8;
  var min_opacity = 0.2

  d3.select(cgm.clust.params.root+' .enr_wait_circle')
    .transition()
    .ease('linear')
    .style('opacity', min_opacity)
    .transition()
    .ease('linear')
    .duration(repeat_time)
    .style('opacity', max_opacity)
    .transition()
    .ease('linear')
    .duration(repeat_time)
    .style('opacity', min_opacity)
    .each("end", animate_wait);
}